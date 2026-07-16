#!/usr/bin/env python3

import argparse
import csv
import re

import pysam


from shared import (
    format_genotype,
    get_alt_depth,
    get_csq_fields,
    get_gene_symbol,
    get_info_value,
    get_sample_depth,
    get_sample_gq,
    has_value,
    load_consequence_order,
    normalize_consequence_term,
    parse_csq_annotations,
    select_consequence,
    select_ensembl_gene_id,
    select_primary_csq,
    value_as_float,
    value_as_text,
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


# Gemini-compatible values consumed by annotate_compound_hets.py.
def gemini_sample_name(sample_name):
    if sample_name in ("0", "-9"):
        return sample_name
    return SAMPLE_NAME_PATTERN.sub("_", str(sample_name))


def get_clinvar_values(record):
    values = []
    for field in ("clinvar_pathogenic", "clinvar_sig", "clinvar_sig_conf"):
        value = get_info_value(record, field)
        if not has_value(value):
            continue
        if isinstance(value, tuple):
            values.extend(str(item) for item in value if has_value(item))
        else:
            values.append(str(value))
    return values


def has_common_clinvar_rescue_value(record):
    return any(value in COMMON_CLINVAR_RESCUE_VALUES for value in get_clinvar_values(record))


def sample_alt_depths(record):
    depths = []
    for sample_name in record.samples:
        alt_depth = get_alt_depth(record.samples[sample_name])
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

    depth = get_sample_depth(sample_data)
    if not has_value(depth):
        return None
    try:
        return int(depth)
    except Exception:
        return depth


def sample_gq_value(sample_data):
    gq = get_sample_gq(sample_data)
    if not has_value(gq):
        return None
    try:
        return float(gq)
    except Exception:
        return gq


def ad_rule(depths, threshold):
    return any(depth >= threshold for depth in depths) or all(depth == -1 for depth in depths)


def gt_type(sample_data):
    """Reproduce Gemini/vcf2db gt_types (0 hom-ref, 1 het, 2 missing, 3 hom-alt).

    compound_hets.py selects candidates on gt_types == 1, so this column decides which
    variants enter the CH analysis and has to agree with the CRE TSVs exactly. Two of
    Gemini's rules are artifacts we reproduce deliberately rather than correct:
      - a hemizygous X/Y call is reported as het, not hom-alt
      - a ref-only half-call is hom-ref when unphased but het when phased, even though
        neither carries an alt read
    """
    gt = sample_data.get("GT")
    if gt is None:
        return 2

    called = [allele for allele in gt if allele is not None]

    if not called:
        # A diploid phased no-call (.|.) is exported as hom-alt; a haploid one is missing.
        # pysam reports phased=True for any single-allele call, hence the ploidy check.
        # annotate_compound_hets.py resets GT == "./." to missing, so this only affects
        # raw-TSV parity, not CH results.
        return 3 if sample_data.phased and len(gt) > 1 else 2

    if len(called) != len(gt):
        # Half-call. An alt allele makes it het. A lone ref allele is hom-ref, except that
        # Gemini calls it het when the genotype is phased and the missing allele is not
        # the leading one (0|. but not .|0) -- despite neither carrying an alt read.
        if any(allele != 0 for allele in called):
            return 1
        return 1 if sample_data.phased and gt[0] is not None else 0

    if all(allele == 0 for allele in called):
        return 0

    if len(called) == 2 and called[0] == called[1]:
        return 3

    return 1


def gt_phase(sample_data):
    return "1" if sample_data.phased else "0"


def record_ps(record, samples):
    value = get_info_value(record, "PS")
    if isinstance(value, tuple):
        parts = list(value)
    elif value is None:
        parts = []
    else:
        parts = str(value).split(",")

    normalized = []
    for idx in range(len(samples)):
        item = parts[idx] if idx < len(parts) else "."
        normalized.append(str(item) if has_value(item) else ".")
    return ",".join(normalized)


def selected_is_high_med(primary, order_map):
    if primary is None:
        return False
    cutoff = order_map.get("IMPACTFUL_CUTOFF")
    if cutoff is None:
        raise SystemExit("IMPACTFUL_CUTOFF not found in impact order file")
    consequence = select_consequence(primary, order_map)
    rank = order_map.get(normalize_consequence_term(consequence), 10**9)
    return rank < cutoff


def common_pathogenic_clinvar(record, faf):
    if faf is None or faf <= MAX_AF:
        return False
    clinvar_status = get_info_value(record, "clinvar_status")
    if not has_value(clinvar_status):
        return False
    if str(clinvar_status) == "no_assertion_criteria_provided":
        return False
    return has_common_clinvar_rescue_value(record)


def output_groups(record, primary, order_map):
    faf = value_as_float(get_info_value(record, "gnomad_af_grpmax"))
    depths = sample_alt_depths(record)
    groups = set()

    is_rare_or_missing = faf is None or faf <= MAX_AF

    if is_rare_or_missing and ad_rule(depths, 3):
        groups.add(HIGH_MED if selected_is_high_med(primary, order_map) else LOW)

    if is_rare_or_missing and get_clinvar_values(record) and ad_rule(depths, 1):
        groups.update([HIGH_MED, LOW])

    if common_pathogenic_clinvar(record, faf):
        groups.update([HIGH_MED, LOW])

    return groups


def row_for_record(record, csq_records, primary, samples, order_map):
    ref = record.ref
    alt = record.alts[0] if record.alts else ""
    variation = select_consequence(primary, order_map) if primary is not None else ""
    gnomad_af_grpmax = value_as_float(get_info_value(record, "gnomad_af_grpmax"))

    row = {
        "Chrom": record.chrom,
        "Pos": record.pos,
        "Variant_id": f"{record.chrom}-{record.pos}-{ref}-{alt}",
        "Ref": ref,
        "Alt": alt,
        "Variation": variation,
        "Depth": value_as_text(get_info_value(record, "DP")),
        "Quality": "" if record.qual is None else record.qual,
        "Gene": get_gene_symbol(primary) if primary is not None else "",
        "Clinvar": ";".join(get_clinvar_values(record)),
        "Ensembl_gene_id": select_ensembl_gene_id(csq_records, primary),
        "Gnomad_af_grpmax": -1 if gnomad_af_grpmax is None else gnomad_af_grpmax,
        "Cadd_score": value_as_text(get_info_value(record, "CADD_phred")),
        "SpliceAI_score": value_as_text(get_info_value(record, "spliceai_score")),
        "promoterAI_score": value_as_text(get_info_value(record, "promoterAI")),
        "PS": record_ps(record, samples),
        "Nucleotide_change_ensembl": primary.get("HGVSc", "") if primary is not None else "",
        "Protein_change_ensembl": primary.get("HGVSp", "") if primary is not None else "",
    }

    for sample in samples:
        sample_data = record.samples[sample]
        sample_column = gemini_sample_name(sample)
        alt_depth = get_alt_depth(sample_data)
        total_depth = sample_total_depth_value(sample_data)
        genotype_quality = sample_gq_value(sample_data)
        row[f"gts.{sample_column}"] = format_genotype(sample_data, ref, record.alts or [])
        row[f"gt_types.{sample_column}"] = gt_type(sample_data)
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
    order_map = load_consequence_order(args.impact_order_file)

    with pysam.VariantFile(args.vcf) as vcf:
        samples = list(vcf.header.samples)
        csq_fields = get_csq_fields(vcf)
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
            for record in vcf:
                if record.alts is None or len(record.alts) == 0:
                    continue
                filters = list(record.filter.keys())
                if filters and "PASS" not in filters:
                    continue
                if record.alts[0] == "*":
                    continue

                csq_records = parse_csq_annotations(record, csq_fields)
                primary = select_primary_csq(csq_records, order_map)
                if primary is None:
                    continue

                groups = output_groups(record, primary, order_map)
                if not groups:
                    continue

                row = row_for_record(record, csq_records, primary, samples, order_map)
                key = row["Variant_id"]
                for group in sorted(groups):
                    if key in seen[group]:
                        continue
                    seen[group].add(key)
                    writers[group].writerow(row)


if __name__ == "__main__":
    main()
