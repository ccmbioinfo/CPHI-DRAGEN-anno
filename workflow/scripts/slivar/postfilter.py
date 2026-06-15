#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict

import pysam

from shared import (
    MISSING,
    as_float,
    consequence_terms,
    info,
    load_impact_order,
    normalize_consequence_term,
    parse_csq_header,
    parse_csq_records,
    present,
)


def stringify(value):
    if value is None:
        return ""
    if isinstance(value, tuple):
        return ",".join(str(v) for v in value if present(v))
    return str(value)


def values_as_list(value):
    if not present(value):
        return []
    if isinstance(value, tuple):
        return [str(v) for v in value if present(v)]
    return [str(value)]


def contains_pathogenic(value):
    if not present(value):
        return False
    if isinstance(value, tuple):
        return any(contains_pathogenic(v) for v in value)
    return "pathogenic" in str(value).lower()


def clinvar_values(record):
    values = []
    for field in ("clinvar_pathogenic", "clinvar_sig", "clinvar_sig_conf"):
        values.extend(values_as_list(info(record, field)))
    return values


def clinvar_text(record):
    return ";".join(clinvar_values(record))


def has_clinvar(record):
    return len(clinvar_values(record)) > 0


def is_pass(record):
    filters = list(record.filter.keys())
    return len(filters) == 0 or filters == ["PASS"] or "PASS" in filters


def is_star_alt(record):
    return record.alts is not None and len(record.alts) > 0 and record.alts[0] == "*"


def variant_key(record):
    return f"{record.chrom}:{record.pos}:{record.ref}:{record.alts[0]}"


def sample_alt_depths(record):
    depths = []
    for sample_name in record.samples:
        sample = record.samples[sample_name]
        if "AD" not in sample:
            depths.append(-1)
            continue
        ad = sample.get("AD")
        if ad is None or len(ad) < 2:
            depths.append(-1)
            continue
        alt_depth = ad[1]
        if alt_depth is None:
            depths.append(-1)
            continue
        try:
            depths.append(int(alt_depth))
        except Exception:
            depths.append(-1)
    return depths


def ad_rule(depths, threshold):
    return any(depth >= threshold for depth in depths) or all(depth == -1 for depth in depths)


def final_report_ad_rule(depths, threshold):
    return any(depth >= threshold for depth in depths)


def base_exclusion_reason(record):
    if record.alts is None or len(record.alts) == 0:
        return False, "missing_alt"
    if not is_pass(record):
        return False, "not_pass"
    if is_star_alt(record):
        return False, "star_alt"
    return True, ""


def is_impactful_non_low(record):
    severity = info(record, "impact_severity")
    if present(severity):
        if isinstance(severity, tuple):
            return any(str(s).upper() != "LOW" for s in severity if present(s))
        return str(severity).upper() != "LOW"

    impactful = info(record, "impactful")
    if impactful is None:
        return False
    if isinstance(impactful, bool):
        return impactful
    if isinstance(impactful, tuple):
        return any(str(v) not in {"0", "False", "false", "", "."} for v in impactful)
    return str(impactful) not in {"0", "False", "false", "", "."}


def pass_rare_impactful(record, max_af):
    ok, reason = base_exclusion_reason(record)
    if not ok:
        return False, reason

    faf = as_float(info(record, "gnomad_fafmax_faf95_max"))
    if faf is not None and faf > max_af:
        return False, "faf_gt_max"
    if not is_impactful_non_low(record):
        return False, "not_impactful_non_low"

    depths = sample_alt_depths(record)
    if not ad_rule(depths, 3):
        return False, "ad_lt_3"

    return (True, "rare_impactful_ad3_missing_faf") if faf is None else (True, "rare_impactful_ad3")


def pass_rare_main(record, max_af):
    ok, reason = base_exclusion_reason(record)
    if not ok:
        return False, reason

    faf = as_float(info(record, "gnomad_fafmax_faf95_max"))
    if faf is not None and faf > max_af:
        return False, "faf_gt_main_max"

    depths = sample_alt_depths(record)
    if not ad_rule(depths, 3):
        return False, "ad_lt_3"

    return (True, "rare_main_ad3_missing_faf") if faf is None else (True, "rare_main_ad3")


def pass_rare_clinvar(record, max_af, allow_missing_faf):
    ok, reason = base_exclusion_reason(record)
    if not ok:
        return False, reason

    faf = as_float(info(record, "gnomad_fafmax_faf95_max"))
    if faf is None:
        if not allow_missing_faf:
            return False, "missing_faf"
    elif faf > max_af:
        return False, "faf_gt_clinvar_max"

    if not has_clinvar(record):
        return False, "no_clinvar"

    depths = sample_alt_depths(record)
    if not ad_rule(depths, 1):
        return False, "ad_lt_1"

    return (True, "rare_clinvar_ad1_missing_faf") if faf is None else (True, "rare_clinvar_ad1")


def pass_common_pathogenic_clinvar(record, common_min_af):
    ok, reason = base_exclusion_reason(record)
    if not ok:
        return False, reason

    faf = as_float(info(record, "gnomad_fafmax_faf95_max"))
    if faf is None:
        return False, "missing_faf"
    if faf < common_min_af:
        return False, "faf_lt_common_min"

    clinvar_status = info(record, "clinvar_status")
    if not present(clinvar_status):
        return False, "missing_clinvar_status"
    if str(clinvar_status) == "no_assertion_criteria_provided":
        return False, "no_assertion_criteria_provided"
    if not contains_pathogenic(clinvar_text(record)):
        return False, "clinvar_not_pathogenic"
    return True, "common_pathogenic_clinvar"


def parse_spliceai_score(record):
    value = info(record, "spliceai_score")
    if not present(value):
        return 0.0

    annotations = []
    if isinstance(value, tuple):
        for item in value:
            annotations.extend(str(item).split(","))
    else:
        annotations = str(value).split(",")

    best_score = 0.0
    for annotation in annotations:
        parts = annotation.split("|")
        if len(parts) < 10:
            continue
        for idx in (2, 3, 4, 5):
            try:
                score = float(parts[idx])
            except Exception:
                continue
            best_score = max(best_score, score)
    return best_score


def parse_promoterai_score(record):
    value = info(record, "promoterAI")
    if not present(value):
        return 0.0

    raw_values = value if isinstance(value, tuple) else str(value).split(",")
    parsed = []
    for raw in raw_values:
        text = str(raw)
        if text in MISSING or text in {"No", "no"}:
            continue
        try:
            parsed.append(abs(float(text)))
        except Exception:
            continue
    return max(parsed) if parsed else 0.0


def top_consequence(csq, order_map):
    terms = consequence_terms(csq)
    if not terms:
        return ""
    return sorted(
        terms,
        key=lambda term: (order_map.get(normalize_consequence_term(term), 999), term),
    )[0]


def csq_gene_symbol(csq):
    for field in ("SYMBOL", "HGNC"):
        value = csq.get(field, "")
        if present(value):
            return value
    gene_field = csq.get("Gene", "")
    if present(gene_field) and not str(gene_field).startswith("ENSG"):
        return gene_field
    return ""


def csq_is_pseudogene(csq):
    biotype = csq.get("BIOTYPE", "")
    return present(biotype) and "pseudogene" in str(biotype).lower()


def csq_is_protein_coding(csq):
    return csq.get("BIOTYPE", "") == "protein_coding"


def csq_is_canonical(csq):
    return csq.get("CANONICAL", "") == "YES"


def csq_mane_select(csq):
    value = csq.get("MANE_SELECT", "")
    return value if present(value) else ""


def csq_mane_plus_clinical(csq):
    value = csq.get("MANE_PLUS_CLINICAL", "")
    return value if present(value) else ""


def csq_mane_rank(csq):
    if csq_mane_select(csq):
        return 0
    if csq_mane_plus_clinical(csq):
        return 1
    return 2


def primary_csq_sort_key(csq, order_map, mode):
    consequence = top_consequence(csq, order_map)
    if csq_is_pseudogene(csq):
        biotype_rank = 2
    elif csq_is_protein_coding(csq):
        biotype_rank = 0
    else:
        biotype_rank = 1
    pseudogene_rank = 1 if csq_is_pseudogene(csq) else 0
    gene_symbol_rank = 0 if present(csq_gene_symbol(csq)) else 1
    canonical_rank = 0 if csq_is_canonical(csq) else 1
    if mode == "coding":
        sort_prefix = (consequence_rank := order_map.get(normalize_consequence_term(consequence), 999), canonical_rank, biotype_rank)
    else:
        sort_prefix = (biotype_rank, consequence_rank := order_map.get(normalize_consequence_term(consequence), 999), gene_symbol_rank, canonical_rank)
    return (
        csq_mane_rank(csq),
        pseudogene_rank,
        *sort_prefix,
        csq.get("Feature", ""),
        csq.get("_csq_index", 0),
    )


def choose_primary_csq(csq_records, order_map, mode):
    if not csq_records:
        return None
    return sorted(csq_records, key=lambda csq: primary_csq_sort_key(csq, order_map, mode))[0]


def all_csq_intergenic(csq_records):
    if not csq_records:
        return False
    saw_any_term = False
    for csq in csq_records:
        terms = consequence_terms(csq)
        if not terms:
            continue
        saw_any_term = True
        if any(term != "intergenic_variant" for term in terms):
            return False
    return saw_any_term


def pass_high_impact_r_style(record, csq_fields, order_map, mode):
    spliceai_score = parse_spliceai_score(record)
    cadd_score = info(record, "CADD_phred")
    promoterai_score = parse_promoterai_score(record)

    cadd_is_missing = not present(cadd_score)
    cadd_numeric = as_float(cadd_score)
    passes_score_gate = (
        spliceai_score >= 0.2
        or cadd_is_missing
        or (cadd_numeric is not None and cadd_numeric >= 10)
        or promoterai_score >= 0.1
    )
    if not passes_score_gate:
        return False, "fails_high_impact_score_gate"

    csq_records = parse_csq_records(record, csq_fields)
    if all_csq_intergenic(csq_records):
        return False, "all_csq_intergenic"

    primary = choose_primary_csq(csq_records, order_map, mode)
    if primary is not None and not present(csq_gene_symbol(primary)):
        return False, "missing_primary_gene_symbol"
    return True, "passes_numeric_score_gate"


def evaluate_record(mode, record, branch, csq_fields, order_map):
    if mode == "coding":
        if branch == "rare_impactful":
            return pass_rare_impactful(record, 0.01)
        if branch == "rare_clinvar":
            return pass_rare_clinvar(record, 0.01, allow_missing_faf=False)
        if branch == "rare_clinvar_allow_missing_faf":
            return pass_rare_clinvar(record, 0.01, allow_missing_faf=True)
        if branch == "common_pathogenic_clinvar":
            return pass_common_pathogenic_clinvar(record, 0.01)
        raise ValueError(f"Unknown coding branch: {branch}")

    rare_main_af = 0.001 if mode == "wgs-high-impact" else 0.01

    if branch == "rare_main":
        keep, reason = pass_rare_main(record, rare_main_af)
    elif branch == "rare_clinvar":
        keep, reason = pass_rare_clinvar(record, 0.01, allow_missing_faf=False)
    elif branch == "rare_clinvar_allow_missing_faf":
        keep, reason = pass_rare_clinvar(record, 0.01, allow_missing_faf=True)
    elif branch == "common_pathogenic_clinvar":
        keep, reason = pass_common_pathogenic_clinvar(record, 0.01)
    else:
        raise ValueError(f"Unknown high-impact branch: {branch}")

    if not keep:
        return False, reason
    if mode == "wgs-high-impact":
        return pass_high_impact_r_style(record, csq_fields, order_map, mode)
    return True, reason


def audit_row(mode, record, branch, keep, reason, csq_fields, order_map):
    depths = sample_alt_depths(record)
    filters = list(record.filter.keys())
    base = {
        "key": variant_key(record),
        "branch": branch,
        "keep": keep,
        "reason": reason,
        "chrom": record.chrom,
        "pos": record.pos,
        "ref": record.ref,
        "alt": record.alts[0] if record.alts else "",
        "filter": ";".join(filters) if filters else "PASS",
        "gnomad_fafmax_faf95_max": "" if as_float(info(record, "gnomad_fafmax_faf95_max")) is None else as_float(info(record, "gnomad_fafmax_faf95_max")),
        "clinvar_status": stringify(info(record, "clinvar_status")),
        "clinvar_text": clinvar_text(record),
        "alt_depths": ",".join(str(depth) for depth in depths),
        "max_alt_depth": max(depths) if depths else "",
        "all_alt_depths_missing": all(depth == -1 for depth in depths),
    }

    if mode == "coding":
        base["impact_severity"] = stringify(info(record, "impact_severity"))
        base["impactful"] = stringify(info(record, "impactful"))
        return base

    csq_records = parse_csq_records(record, csq_fields)
    primary = choose_primary_csq(csq_records, order_map, mode)
    base["cadd_phred"] = stringify(info(record, "CADD_phred"))
    base["spliceai_max_score"] = parse_spliceai_score(record)
    base["promoterAI_abs_max"] = parse_promoterai_score(record)
    base["primary_gene_symbol"] = csq_gene_symbol(primary) if primary is not None else ""
    base["primary_consequence"] = top_consequence(primary, order_map) if primary is not None else ""
    return base


def write_outputs(mode, out_prefix, kept_keys, audit_rows, kept_reasons_by_key, branch_reason_counts, duplicate_kept_branch_records, input_counts, kept_input_counts, dropped_input_counts):
    if mode == "coding":
        suffix = "post_gemini_filter"
    elif mode == "wgs-high-impact":
        suffix = "post_high_impact_filter"
    else:
        suffix = "post_wgs_filter"

    with open(f"{out_prefix}.{suffix}.keys.txt", "w") as handle:
        for record_key in sorted(kept_keys):
            handle.write(f"{record_key}\n")

    fields = [
        "key",
        "branch",
        "keep",
        "reason",
        "chrom",
        "pos",
        "ref",
        "alt",
        "filter",
        "gnomad_fafmax_faf95_max",
        "clinvar_status",
        "clinvar_text",
        "alt_depths",
        "max_alt_depth",
        "all_alt_depths_missing",
    ]
    if mode == "coding":
        fields.extend(["impact_severity", "impactful"])
    else:
        fields.extend(["cadd_phred", "spliceai_max_score", "promoterAI_abs_max", "primary_gene_symbol", "primary_consequence"])

    with open(f"{out_prefix}.{suffix}.audit.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(audit_rows)

    with open(f"{out_prefix}.{suffix}.reasons.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["key", "reasons"])
        for record_key in sorted(kept_reasons_by_key):
            writer.writerow([record_key, ",".join(sorted(set(kept_reasons_by_key[record_key])) )])

    with open(f"{out_prefix}.{suffix}.summary.txt", "w") as handle:
        handle.write(f"selected_unique_keys\t{len(kept_keys)}\n")
        handle.write(f"duplicate_kept_branch_records\t{duplicate_kept_branch_records}\n")
        handle.write("\ninput_records_by_branch\n")
        for branch in sorted(input_counts):
            handle.write(f"{branch}\t{input_counts[branch]}\n")
        handle.write("\nkept_records_by_branch\n")
        for branch in sorted(kept_input_counts):
            handle.write(f"{branch}\t{kept_input_counts[branch]}\n")
        handle.write("\ndropped_records_by_branch\n")
        for branch in sorted(dropped_input_counts):
            handle.write(f"{branch}\t{dropped_input_counts[branch]}\n")
        handle.write("\nkept_reason_counts\n")
        for reason in sorted(branch_reason_counts):
            handle.write(f"{reason}\t{branch_reason_counts[reason]}\n")


def parse_args():
    parser = argparse.ArgumentParser(description="Unified slivar post-filtering for coding, wgs, and wgs-high-impact.")
    parser.add_argument("--mode", required=True, choices=["coding", "wgs", "wgs-high-impact"])
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--impact-order-file", required=True)
    parser.add_argument("--rare-impactful-vcf")
    parser.add_argument("--rare-main-vcf")
    parser.add_argument("--rare-clinvar-vcf", required=True)
    parser.add_argument("--common-pathogenic-clinvar-vcf", required=True)
    parser.add_argument("--rare-clinvar-allow-missing-faf-vcf", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    order_map = load_impact_order(args.impact_order_file)

    if args.mode == "coding":
        branch_inputs = [
            ("rare_impactful", args.rare_impactful_vcf),
            ("rare_clinvar", args.rare_clinvar_vcf),
            ("common_pathogenic_clinvar", args.common_pathogenic_clinvar_vcf),
            ("rare_clinvar_allow_missing_faf", args.rare_clinvar_allow_missing_faf_vcf),
        ]
        out_vcf_path = f"{args.out_prefix}.post_gemini_filter.vcf"
    elif args.mode == "wgs-high-impact":
        branch_inputs = [
            ("rare_main", args.rare_main_vcf),
            ("rare_clinvar", args.rare_clinvar_vcf),
            ("common_pathogenic_clinvar", args.common_pathogenic_clinvar_vcf),
            ("rare_clinvar_allow_missing_faf", args.rare_clinvar_allow_missing_faf_vcf),
        ]
        out_vcf_path = f"{args.out_prefix}.post_high_impact_filter.vcf"
    else:
        branch_inputs = [
            ("rare_main", args.rare_main_vcf),
            ("rare_clinvar", args.rare_clinvar_vcf),
            ("common_pathogenic_clinvar", args.common_pathogenic_clinvar_vcf),
            ("rare_clinvar_allow_missing_faf", args.rare_clinvar_allow_missing_faf_vcf),
        ]
        out_vcf_path = f"{args.out_prefix}.post_wgs_filter.vcf"

    with pysam.VariantFile(branch_inputs[0][1]) as first_vcf:
        csq_fields = parse_csq_header(first_vcf)
        if args.mode in {"wgs", "wgs-high-impact"} and not csq_fields:
            raise SystemExit("VCF header does not contain a usable INFO/CSQ definition")
        out_vcf = pysam.VariantFile(out_vcf_path, "w", header=first_vcf.header)

        seen_kept_keys = set()
        kept_keys = []
        audit_rows = []
        kept_reasons_by_key = defaultdict(list)
        branch_reason_counts = defaultdict(int)
        duplicate_kept_branch_records = 0
        input_counts = defaultdict(int)
        kept_input_counts = defaultdict(int)
        dropped_input_counts = defaultdict(int)

        for branch, path in branch_inputs:
            with pysam.VariantFile(path) as branch_vcf:
                for record in branch_vcf:
                    if record.alts is None or len(record.alts) == 0:
                        continue

                    record_key = variant_key(record)
                    input_counts[branch] += 1
                    keep, reason = evaluate_record(args.mode, record, branch, csq_fields, order_map)
                    if keep:
                        depths = sample_alt_depths(record)
                        if not final_report_ad_rule(depths, 3):
                            keep, reason = False, "final_report_ad_lt_3"
                    audit_rows.append(audit_row(args.mode, record, branch, keep, reason, csq_fields, order_map))

                    if not keep:
                        dropped_input_counts[branch] += 1
                        continue

                    kept_input_counts[branch] += 1
                    branch_reason_counts[reason] += 1
                    kept_reasons_by_key[record_key].append(reason)

                    if record_key in seen_kept_keys:
                        duplicate_kept_branch_records += 1
                        continue

                    seen_kept_keys.add(record_key)
                    kept_keys.append(record_key)
                    out_vcf.write(record)

        out_vcf.close()

    write_outputs(
        args.mode,
        args.out_prefix,
        kept_keys,
        audit_rows,
        kept_reasons_by_key,
        branch_reason_counts,
        duplicate_kept_branch_records,
        input_counts,
        kept_input_counts,
        dropped_input_counts,
    )


if __name__ == "__main__":
    main()
