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


# Get the highest SpliceAI score.
def parse_spliceai_score(record):
    value = get_info_value(record, "spliceai_score")
    if not has_value(value):
        return 0.0

    # Roughly follows the SpliceAI parsing in workflow/scripts/cre/cre.vcf2db.R.
    # Values look like Allele|Gene|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|...
    # Repeated INFO values might be a tuple. Each value can also
    # still hold comma-separated annotations, so flatten both forms first.
    annotations = []
    if isinstance(value, tuple):
        for item in value:
            annotations.extend(str(item).split(","))
    else:
        annotations = str(value).split(",")

    # Use the highest delta score across every SpliceAI annotation.
    best_score = 0.0
    for annotation in annotations:
        parts = annotation.split("|")
        # Short or malformed annotations cannot contain all four DS values.
        if len(parts) < 10:
            continue
        # SpliceAI DS_AG/DS_AL/DS_DG/DS_DL occupy fields 2-5 in this annotation.
        for idx in (2, 3, 4, 5):
            try:
                score = float(parts[idx])
            except Exception:
                continue
            best_score = max(best_score, score)
    return best_score


# Get the highest absolute promoterAI score.
def parse_promoterai_score(record):
    value = get_info_value(record, "promoterAI")
    if not has_value(value):
        return 0.0

    # Roughly follows the promoterAI parsing in workflow/scripts/cre/cre.vcf2db.R.
    # Values are numeric score strings, for example "0.22" or "-0.14,0.03".
    # Repeated INFO values might be a tuple. Each value can also hold
    # comma-separated scores, so normalize to an iterable first.
    if isinstance(value, tuple):
        raw_values = value
    else:
        raw_values = str(value).split(",")
    parsed = []
    for raw in raw_values:
        text = str(raw)
        # CRE treated missing values and "No" promoterAI hits as zero.
        if text in MISSING or text in {"No", "no"}:
            continue
        try:
            parsed.append(abs(float(text)))
        except Exception:
            continue
    if parsed:
        return max(parsed)
    return 0.0


# True only when at least one parsed consequence exists and all are intergenic.
def all_csq_intergenic(csq_records):
    if not csq_records:
        return False
    all_intergenic = False
    for csq in csq_records:
        terms = get_consequence_terms(csq)
        if not terms:
            continue
        if any(term != "intergenic_variant" for term in terms):
            return False
        all_intergenic = True
    return all_intergenic

# For a given CSQ, return sort keys used only by the high-impact inclusion gate.
def high_impact_csq_sort_key(csq, consequence_order):
    biotype = str(csq.get("BIOTYPE", "")).lower()

    if has_value(csq.get("MANE_SELECT", "")):
        mane_rank = 0
    elif has_value(csq.get("MANE_PLUS_CLINICAL", "")):
        mane_rank = 1
    else:
        mane_rank = 2

    if "pseudogene" in biotype:
        pseudogene_rank = 1
    else:
        pseudogene_rank = 0

    # Favors ordinary protein-coding annotations and de-prioritizes
    # pseudogenes before checking whether the selected CSQ has a gene symbol.
    if biotype == "protein_coding":
        biotype_rank = 0
    elif "pseudogene" in biotype:
        biotype_rank = 2
    else:
        biotype_rank = 1

    if has_value(get_gene_symbol(csq)):
        gene_rank = 0
    else:
        gene_rank = 1

    if csq.get("CANONICAL", "") == "YES":
        canonical_rank = 0
    else:
        canonical_rank = 1

    # min() compares tuple entries from left to right, so lower ranks win in
    # this order: MANE, non-pseudogene, biotype, consequence, gene, canonical.
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


# Retain a high-impact variant only when it passes the noncoding score gate,
# is not annotated as entirely intergenic, and the most impactful CSQ chosen
# for this inclusion gate has an associated gene symbol. This CSQ choice is
# separate from the final report transcript selection in build_report.py.
def passes_high_impact_gate(record, csq_fields, consequence_order):
    spliceai_score = parse_spliceai_score(record)
    cadd_score = get_info_value(record, "CADD_phred")
    promoterai_score = parse_promoterai_score(record)

    # Match the CRE high-impact score gate: SpliceAI >= 0.2, missing CADD,
    # CADD >= 10, or absolute promoterAI >= 0.1.
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

    # Parse CSQ into field-name dictionaries for consequence and gene checks.
    csq_records = parse_csq_annotations(record, csq_fields)
    # Drop records only when every annotated consequence is explicitly intergenic.
    if all_csq_intergenic(csq_records):
        return False

    # This primary CSQ is only for the high-impact inclusion test; build_report.py
    # performs the final report-display CSQ selection.
    primary = None
    if csq_records:
        # Choose the CSQ that sorts first by the high-impact gate order above.
        primary = min(
            csq_records,
            key=lambda csq: high_impact_csq_sort_key(csq, consequence_order),
        )
    # Final check to ensure there is a gene symbol for the chosen primary CSQ.
    if primary is None or not has_value(get_gene_symbol(primary)):
        return False
    return True


def parse_args():
    parser = argparse.ArgumentParser(description="Post-filter and merge slivar report branches.")
    parser.add_argument(
        "--mode",
        required=True,
        choices=["coding", "wgs-high-impact", "denovo", "panel", "panel-flank"],
    )
    parser.add_argument("--out-vcf", required=True)
    parser.add_argument("--impact-order-file", required=True)
    parser.add_argument("--rare-main-vcf", required=True)
    parser.add_argument("--rare-clinvar-vcf", required=True)
    parser.add_argument("--common-pathogenic-clinvar-vcf", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    is_high_impact = args.mode == "wgs-high-impact"
    consequence_order = None
    if is_high_impact:
        # The high-impact gene exists gate ranks CSQ terms.
        consequence_order = load_consequence_order(args.impact_order_file)

    branch_inputs = [
        args.rare_main_vcf,
        args.rare_clinvar_vcf,
        args.common_pathogenic_clinvar_vcf,
    ]

    # slivar expr has already applied each branch's selection criteria. This step only
    # merges and deduplicates those branches, plus the additional high-impact gate.
    with pysam.VariantFile(branch_inputs[0]) as first_vcf:
        if is_high_impact:
            # High-impact filtering parses VEP CSQ annotations by field name.
            csq_fields = get_csq_fields(first_vcf)
        else:
            csq_fields = []
        if is_high_impact and not csq_fields:
            raise SystemExit("VCF header does not contain a usable INFO/CSQ definition")
        # The output VCF uses the rare-main branch header.
        out_vcf = pysam.VariantFile(args.out_vcf, "w", header=first_vcf.header)

        # Walk the branch VCFs in order and write each new variant key once.
        # Could have used bcftools merge if not for high-impact require some custom gates.
        seen_kept_keys = set()

        for path in branch_inputs:
            with pysam.VariantFile(path) as branch_vcf:
                for record in branch_vcf:
                    # Branch inputs are decomposed; use the first ALT as the key allele.
                    if record.alts:
                        alt = record.alts[0]
                    else:
                        alt = ""
                    # Compare variants based on key.
                    record_key = (
                        record.chrom,
                        record.pos,
                        record.ref,
                        alt,
                    )
                    # Only high-impact mode needs the extra score/intergenic/gene-symbol gate.
                    if is_high_impact and not passes_high_impact_gate(record, csq_fields, consequence_order):
                        continue

                    if record_key in seen_kept_keys:
                        continue

                    seen_kept_keys.add(record_key)
                    out_vcf.write(record)

        out_vcf.close()


if __name__ == "__main__":
    main()
