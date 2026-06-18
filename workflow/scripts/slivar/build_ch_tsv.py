#!/usr/bin/env python3

import argparse
import csv
import re

import pysam

from shared import (
    as_float,
    as_text,
    choose_primary_csq,
    consequence_display,
    gene_id_with_fallback_for_csq,
    gene_symbol,
    gt_string,
    info,
    load_impact_order,
    normalize_consequence_term,
    parse_csq_header,
    parse_csq_records,
    present,
    sample_alt_depth_value,
    sample_depth,
    sample_gq,
)


MAX_AF = 0.01
HIGH_MED = "HIGH-MED"
LOW = "LOW"
COMMON_CLINVAR_RESCUE_VALUES = {
    "Pathogenic",
    "Likely_pathogenic",
    "Conflicting_interpretations_of_pathogenicity",
}
SAMPLE_NAME_PATTERN = re.compile(r"[-\s\\]")


def gemini_sample_name(sample_name):
    if sample_name in ("0", "-9"):
        return sample_name
    return SAMPLE_NAME_PATTERN.sub("_", str(sample_name))


def values_as_list(value):
    if not present(value):
        return []
    if isinstance(value, tuple):
        return [str(v) for v in value if present(v)]
    return [str(value)]


def clinvar_values(record):
    values = []
    for field in ("clinvar_pathogenic", "clinvar_sig", "clinvar_sig_conf"):
        values.extend(values_as_list(info(record, field)))
    return values


def clinvar_text(record):
    return ";".join(clinvar_values(record))


def has_clinvar(record):
    return len(clinvar_values(record)) > 0


def has_common_clinvar_rescue_value(record):
    return any(value in COMMON_CLINVAR_RESCUE_VALUES for value in clinvar_values(record))


def is_pass(record):
    filters = list(record.filter.keys())
    return len(filters) == 0 or filters == ["PASS"] or "PASS" in filters


def is_star_alt(record):
    return record.alts is not None and len(record.alts) > 0 and record.alts[0] == "*"


def sample_alt_depths(record):
    depths = []
    for sample_name in record.samples:
        alt_depth = sample_alt_depth_value(record.samples[sample_name])
        depths.append(-1 if alt_depth is None else alt_depth)
    return depths


def sample_total_depth_value(sample_data):
    ad = sample_data.get("AD")
    if ad is not None:
        if isinstance(ad, str):
            depths = ad.split(",")
        else:
            try:
                depths = list(ad)
            except Exception:
                depths = []

        parsed = []
        for depth in depths:
            if depth is None or str(depth) in {"", ".", "None", "NA"}:
                continue
            try:
                parsed.append(int(depth))
            except Exception:
                continue
        if parsed:
            return sum(parsed)

    depth = sample_depth(sample_data)
    if not present(depth):
        return None
    try:
        return int(depth)
    except Exception:
        return depth


def sample_gq_value(sample_data):
    gq = sample_gq(sample_data)
    if not present(gq):
        return None
    try:
        return float(gq)
    except Exception:
        return gq


def ad_rule(depths, threshold):
    return any(depth >= threshold for depth in depths) or all(depth == -1 for depth in depths)


def gt_type(sample_data, chrom=""):
    gt = sample_data.get("GT")
    if gt is None or all(allele is None for allele in gt):
        return 2

    alt_depth = sample_alt_depth_value(sample_data) or 0
    if alt_depth == 0:
        return 0

    called = [allele for allele in gt if allele is not None]
    if not called:
        return 2

    if all(allele == 0 for allele in called):
        return 0

    if len(called) == 1 and called[0] != 0:
        chrom_text = str(chrom)
        if "X" in chrom_text or "Y" in chrom_text:
            return 3

    if len(called) == 2 and all(allele == called[0] for allele in called) and called[0] != 0:
        return 3

    return 1


def gt_phase(sample_data):
    return "1" if sample_data.phased else "0"


def record_ps(record, samples):
    value = info(record, "PS")
    if isinstance(value, tuple):
        parts = list(value)
    elif value is None:
        parts = []
    else:
        parts = str(value).split(",")

    normalized = []
    for idx in range(len(samples)):
        item = parts[idx] if idx < len(parts) else "."
        normalized.append(str(item) if present(item) else ".")
    return ",".join(normalized)


def selected_is_high_med(primary, order_map):
    if primary is None:
        return False
    cutoff = order_map.get("IMPACTFUL_CUTOFF")
    if cutoff is None:
        raise SystemExit("IMPACTFUL_CUTOFF not found in impact order file")
    consequence = consequence_display(primary, order_map)
    rank = order_map.get(normalize_consequence_term(consequence), 10**9)
    return rank < cutoff


def common_pathogenic_clinvar(record, faf):
    if faf is None or faf <= MAX_AF:
        return False
    clinvar_status = info(record, "clinvar_status")
    if not present(clinvar_status):
        return False
    if str(clinvar_status) == "no_assertion_criteria_provided":
        return False
    return has_common_clinvar_rescue_value(record)


def output_groups(record, primary, order_map):
    faf = as_float(info(record, "gnomad_af_grpmax"))
    depths = sample_alt_depths(record)
    groups = set()

    is_rare_or_missing = faf is None or faf <= MAX_AF

    if is_rare_or_missing and ad_rule(depths, 3):
        groups.add(HIGH_MED if selected_is_high_med(primary, order_map) else LOW)

    if is_rare_or_missing and has_clinvar(record) and ad_rule(depths, 1):
        groups.update([HIGH_MED, LOW])

    if common_pathogenic_clinvar(record, faf):
        groups.update([HIGH_MED, LOW])

    return groups


def row_for_record(record, csq_records, primary, samples, variant_index, order_map):
    ref = record.ref
    alt = record.alts[0] if record.alts else ""
    variation = consequence_display(primary, order_map) if primary is not None else ""
    gnomad_af_grpmax = as_float(info(record, "gnomad_af_grpmax"))

    row = {
        "Chrom": record.chrom,
        "Pos": record.pos,
        "Variant_id": f"{record.chrom}-{record.pos}-{ref}-{alt}",
        "Ref": ref,
        "Alt": alt,
        "Variation": variation,
        "Depth": as_text(info(record, "DP")),
        "Quality": "" if record.qual is None else record.qual,
        "Gene": gene_symbol(primary) if primary is not None else "",
        "Clinvar": clinvar_text(record),
        "Ensembl_gene_id": gene_id_with_fallback_for_csq(csq_records, primary),
        "Gnomad_af_grpmax": -1 if gnomad_af_grpmax is None else gnomad_af_grpmax,
        "Cadd_score": as_text(info(record, "CADD_phred")),
        "SpliceAI_score": as_text(info(record, "spliceai_score")),
        "promoterAI_score": as_text(info(record, "promoterAI")),
        "PS": record_ps(record, samples),
        "Nucleotide_change_ensembl": primary.get("HGVSc", "") if primary is not None else "",
        "Protein_change_ensembl": primary.get("HGVSp", "") if primary is not None else "",
    }

    for sample in samples:
        sample_data = record.samples[sample]
        sample_column = gemini_sample_name(sample)
        alt_depth = sample_alt_depth_value(sample_data)
        total_depth = sample_total_depth_value(sample_data)
        genotype_quality = sample_gq_value(sample_data)
        row[f"gts.{sample_column}"] = gt_string(sample_data, ref, record.alts or [])
        row[f"gt_types.{sample_column}"] = gt_type(sample_data, record.chrom)
        row[f"gt_phases.{sample_column}"] = gt_phase(sample_data)
        row[f"gt_alt_depths.{sample_column}"] = -1 if alt_depth is None else alt_depth
        row[f"gt_depths.{sample_column}"] = -1 if total_depth is None else total_depth
        row[f"gt_quals.{sample_column}"] = -1 if genotype_quality is None else genotype_quality

    return row


def fieldnames(samples):
    fields = [
        "Chrom",
        "Pos",
        "Variant_id",
        "Ref",
        "Alt",
        "Variation",
        "Depth",
        "Quality",
        "Gene",
        "Clinvar",
        "Ensembl_gene_id",
        "Gnomad_af_grpmax",
        "Cadd_score",
        "SpliceAI_score",
        "promoterAI_score",
        "PS",
        "Nucleotide_change_ensembl",
        "Protein_change_ensembl",
    ]
    for sample in samples:
        sample_column = gemini_sample_name(sample)
        fields.extend(
            [
                f"gts.{sample_column}",
                f"gt_types.{sample_column}",
                f"gt_phases.{sample_column}",
                f"gt_alt_depths.{sample_column}",
                f"gt_depths.{sample_column}",
                f"gt_quals.{sample_column}",
            ]
        )
    return fields


def parse_args():
    parser = argparse.ArgumentParser(description="Build slivar-derived sequence variant TSVs for compound-het annotation.")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--impact-order-file", required=True)
    parser.add_argument("--high-med-out", required=True)
    parser.add_argument("--low-out", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    order_map = load_impact_order(args.impact_order_file)

    with pysam.VariantFile(args.vcf) as vcf:
        samples = list(vcf.header.samples)
        csq_fields = parse_csq_header(vcf)
        if not csq_fields:
            raise SystemExit("VCF header does not contain a usable INFO/CSQ definition")

        fields = fieldnames(samples)
        with open(args.high_med_out, "w", newline="") as high_handle, open(args.low_out, "w", newline="") as low_handle:
            writers = {
                HIGH_MED: csv.DictWriter(high_handle, delimiter="\t", fieldnames=fields),
                LOW: csv.DictWriter(low_handle, delimiter="\t", fieldnames=fields),
            }
            for writer in writers.values():
                writer.writeheader()

            seen = {HIGH_MED: set(), LOW: set()}
            for variant_index, record in enumerate(vcf, start=1):
                if record.alts is None or len(record.alts) == 0:
                    continue
                if not is_pass(record) or is_star_alt(record):
                    continue

                csq_records = parse_csq_records(record, csq_fields)
                primary = choose_primary_csq(csq_records, order_map, "coding")
                if primary is None:
                    continue

                groups = output_groups(record, primary, order_map)
                if not groups:
                    continue

                row = row_for_record(record, csq_records, primary, samples, variant_index, order_map)
                key = row["Variant_id"]
                for group in sorted(groups):
                    if key in seen[group]:
                        continue
                    seen[group].add(key)
                    writers[group].writerow(row)


if __name__ == "__main__":
    main()
