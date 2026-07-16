#!/usr/bin/env python3

import argparse

import pysam

from shared import (
    MISSING,
    get_consequence_terms,
    get_csq_fields,
    get_gene_symbol,
    get_info_value,
    has_value,
    load_consequence_order,
    parse_csq_annotations,
    rank_consequence,
    value_as_float,
)


# Extra high-impact filtering that is not part of the slivar expr branches.
def parse_spliceai_score(record):
    value = get_info_value(record, "spliceai_score")
    if not has_value(value):
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
    value = get_info_value(record, "promoterAI")
    if not has_value(value):
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


def all_csq_intergenic(csq_records):
    if not csq_records:
        return False
    saw_any_term = False
    for csq in csq_records:
        terms = get_consequence_terms(csq)
        if not terms:
            continue
        saw_any_term = True
        if any(term != "intergenic_variant" for term in terms):
            return False
    return saw_any_term


def high_impact_csq_sort_key(csq, consequence_order):
    """Legacy CSQ ordering used only by the high-impact inclusion gate."""
    biotype = str(csq.get("BIOTYPE", "")).lower()

    if has_value(csq.get("MANE_SELECT", "")):
        mane_rank = 0
    elif has_value(csq.get("MANE_PLUS_CLINICAL", "")):
        mane_rank = 1
    else:
        mane_rank = 2

    pseudogene_rank = 1 if "pseudogene" in biotype else 0

    if biotype == "protein_coding":
        biotype_rank = 0
    elif "pseudogene" in biotype:
        biotype_rank = 2
    else:
        biotype_rank = 1

    gene_rank = 0 if has_value(get_gene_symbol(csq)) else 1
    canonical_rank = 0 if csq.get("CANONICAL", "") == "YES" else 1

    return (
        mane_rank,
        pseudogene_rank,
        biotype_rank,
        rank_consequence(csq, consequence_order),
        gene_rank,
        canonical_rank,
        csq.get("Feature", ""),
        csq.get("_csq_index", 0),
    )


def passes_high_impact_gate(record, csq_fields, consequence_order):
    spliceai_score = parse_spliceai_score(record)
    cadd_score = get_info_value(record, "CADD_phred")
    promoterai_score = parse_promoterai_score(record)

    cadd_is_missing = not has_value(cadd_score)
    cadd_numeric = value_as_float(cadd_score)
    passes_score_gate = (
        spliceai_score >= 0.2
        or cadd_is_missing
        or (cadd_numeric is not None and cadd_numeric >= 10)
        or promoterai_score >= 0.1
    )
    if not passes_score_gate:
        return False

    csq_records = parse_csq_annotations(record, csq_fields)
    if all_csq_intergenic(csq_records):
        return False

    primary = min(
        csq_records,
        key=lambda csq: high_impact_csq_sort_key(csq, consequence_order),
    ) if csq_records else None
    if primary is None or not has_value(get_gene_symbol(primary)):
        return False
    return True


def parse_args():
    parser = argparse.ArgumentParser(description="Post-filter and merge slivar report branches.")
    parser.add_argument("--mode", required=True, choices=["coding", "wgs", "wgs-high-impact"])
    parser.add_argument("--out-vcf", required=True)
    parser.add_argument("--impact-order-file", required=True)
    parser.add_argument("--rare-main-vcf", required=True)
    parser.add_argument("--rare-clinvar-vcf", required=True)
    parser.add_argument("--common-pathogenic-clinvar-vcf", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    is_high_impact = args.mode == "wgs-high-impact"
    consequence_order = load_consequence_order(args.impact_order_file) if is_high_impact else None

    branch_inputs = [
        args.rare_main_vcf,
        args.rare_clinvar_vcf,
        args.common_pathogenic_clinvar_vcf,
    ]

    # slivar expr has already applied each branch's selection criteria. This step only
    # merges and deduplicates those branches, plus the additional high-impact gate.
    with pysam.VariantFile(branch_inputs[0]) as first_vcf:
        csq_fields = get_csq_fields(first_vcf) if is_high_impact else []
        if is_high_impact and not csq_fields:
            raise SystemExit("VCF header does not contain a usable INFO/CSQ definition")
        out_vcf = pysam.VariantFile(args.out_vcf, "w", header=first_vcf.header)

        seen_kept_keys = set()

        for path in branch_inputs:
            with pysam.VariantFile(path) as branch_vcf:
                for record in branch_vcf:
                    record_key = (
                        record.chrom,
                        record.pos,
                        record.ref,
                        record.alts[0] if record.alts else "",
                    )
                    if is_high_impact and not passes_high_impact_gate(record, csq_fields, consequence_order):
                        continue

                    if record_key in seen_kept_keys:
                        continue

                    seen_kept_keys.add(record_key)
                    out_vcf.write(record)

        out_vcf.close()


if __name__ == "__main__":
    main()
