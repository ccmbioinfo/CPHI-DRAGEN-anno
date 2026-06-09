#!/usr/bin/env python3

import argparse
import csv
from collections import Counter, defaultdict

import pysam

from shared import (
    MISSING,
    as_float,
    as_text,
    clean_slash,
    consequence_terms,
    gene_symbol,
    gt_string,
    index_constraint,
    index_first,
    index_many,
    info,
    load_impact_order,
    load_table,
    normalize_consequence_term,
    normalized_consequence_terms,
    parse_csq_header,
    parse_csq_records,
    present,
    sample_alt_depth,
    sample_depth,
    sample_gq,
    uniq_join,
    variant_key,
    zygosity,
)


IMPACT_RANK = {
    "HIGH": 0,
    "MODERATE": 1,
    "LOW": 2,
    "MODIFIER": 3,
}


def consequence_rank(csq, order_map):
    ranks = [order_map[term] for term in normalized_consequence_terms(csq) if term in order_map]
    return min(ranks) if ranks else 10**9


def ensembl_gene_id(csq):
    gene = csq.get("Gene", "")
    return gene if gene.startswith("ENSG") else ""


def best_ensembl_gene_id_for_symbol(csq_records, symbol):
    candidates = []
    for csq in csq_records:
        if gene_symbol(csq) != symbol:
            continue
        gene = ensembl_gene_id(csq)
        if gene:
            candidates.append(gene)
    return candidates[0] if candidates else ""


def parse_vep_score(value):
    if not present(value):
        return ""
    text = str(value)
    if "(" in text and ")" in text:
        return text.rsplit("(", 1)[1].split(")", 1)[0]
    return text


def build_info_item(csq):
    gene = gene_symbol(csq)
    exon = csq.get("EXON", "") or "NA"
    if not present(exon):
        exon = "NA"
    return f"{gene}:exon{exon}:{csq.get('HGVSc', '')}:{csq.get('HGVSp', '')}".replace("%3D", "=")


def build_info_all(csq_records):
    return uniq_join(build_info_item(csq) for csq in csq_records)


def refseq_change_all(csq_records):
    return uniq_join(refseq_change_for_csq(csq) for csq in csq_records if present(refseq_change_for_csq(csq))) or "NA"


def canonical_rank(csq):
    return 0 if csq.get("CANONICAL", "") == "YES" else 1


def mane_value(csq, field):
    value = csq.get(field, "")
    return value if present(value) else ""


def mane_select_value(csq):
    return mane_value(csq, "MANE_SELECT")


def mane_plus_clinical_value(csq):
    return mane_value(csq, "MANE_PLUS_CLINICAL")


def mane_rank(csq):
    if mane_select_value(csq):
        return 0
    if mane_plus_clinical_value(csq):
        return 1
    return 2


def biotype_rank(csq, mode):
    biotype = csq.get("BIOTYPE", "")
    if biotype == "protein_coding":
        return 0
    if mode != "coding" and "pseudogene" in biotype:
        return 2
    return 1


def gene_symbol_rank(csq):
    return 0 if present(gene_symbol(csq)) else 1


def csq_sort_key(csq, order_map, mode):
    if mode == "coding":
        return (
            mane_rank(csq),
            consequence_rank(csq, order_map),
            IMPACT_RANK.get(csq.get("IMPACT", ""), 99),
            canonical_rank(csq),
            biotype_rank(csq, mode),
            csq.get("Feature", ""),
            csq.get("_csq_index", ""),
        )
    return (
        mane_rank(csq),
        biotype_rank(csq, mode),
        consequence_rank(csq, order_map),
        gene_symbol_rank(csq),
        canonical_rank(csq),
        csq.get("Feature", ""),
        csq.get("_csq_index", ""),
    )


def choose_primary_csq(csq_records, order_map, mode):
    if not csq_records:
        return None
    return sorted(csq_records, key=lambda csq: csq_sort_key(csq, order_map, mode))[0]


def parse_spliceai(value):
    if not present(value):
        return "NA|NA|NA", "0"
    annotations = []
    if isinstance(value, tuple):
        for item in value:
            annotations.extend(str(item).split(","))
    else:
        annotations = str(value).split(",")
    best_gene = "NA"
    best_impact = "NA"
    best_score = 0.0
    best_pos = "NA"
    fields = [("acceptor_gain", 2, 6), ("acceptor_loss", 3, 7), ("donor_gain", 4, 8), ("donor_loss", 5, 9)]
    for annotation in annotations:
        parts = annotation.split("|")
        if len(parts) < 10:
            continue
        gene = parts[1]
        for impact_name, score_idx, pos_idx in fields:
            try:
                score = float(parts[score_idx])
            except Exception:
                continue
            if score > best_score:
                best_score = score
                best_gene = gene
                best_impact = impact_name
                best_pos = parts[pos_idx]
    if best_score == 0:
        return "NA|NA|NA", "0"
    return f"{best_gene}|{best_impact}|{best_pos}", str(best_score)


def parse_promoterai(value):
    if not present(value):
        return "0"
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
    return str(max(parsed)) if parsed else "0"


def noncoding_pred_fraction(cadd, ncer, remm, linsight):
    def float_or_none(value):
        if not present(value):
            return None
        try:
            return float(str(value))
        except Exception:
            return None

    preds = []
    values = [
        (float_or_none(cadd), 10),
        (float_or_none(ncer), 95.95),
        (float_or_none(remm), 0.9585),
        (float_or_none(linsight), 0.9828),
    ]
    for value, threshold in values:
        if value is not None:
            preds.append(value > threshold)
    return "0//0" if not preds else f"{sum(preds)}//{len(preds)}"


def hyperlink_ucsc(position):
    return '=HYPERLINK("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hgt.out3=10x&position=' + position + '","UCSC_link")'


def hyperlink_gnomad(chrom, pos, ref, alt):
    variant_id = f"{chrom}-{pos}-{ref}-{alt}"
    return '=HYPERLINK("http://gnomad.broadinstitute.org/variant/' + variant_id + '?dataset=gnomad_r4","GNOMAD_link")'


def group_records_by_variant(vcf):
    grouped = defaultdict(list)
    for record in vcf:
        if record.alts:
            grouped[variant_key(record)].append(record)
    return grouped


def first_info(records, field, default=""):
    for record in records:
        value = info(record, field)
        if present(value):
            return value
    return default


def merged_clinvar_text(records):
    values = []
    for record in records:
        for field in ("clinvar_pathogenic", "clinvar_sig", "clinvar_sig_conf"):
            value = info(record, field)
            if not present(value):
                continue
            if isinstance(value, tuple):
                values.extend(str(v) for v in value if present(v))
            else:
                values.append(str(value))
    return uniq_join(values, sep=";")


def merged_csq_records(records, csq_fields):
    all_csq = []
    next_index = 1
    for record in records:
        parsed = parse_csq_records(record, csq_fields, start_index=next_index)
        all_csq.extend(parsed)
        next_index += len(parsed)
    return all_csq


def join_gene_description(ensembl_gene_id_value, gene_descriptions_by_ensg):
    return gene_descriptions_by_ensg.get(ensembl_gene_id_value, {}).get("Gene_description", "")


def join_omim(gene, omim_by_gene):
    matches = omim_by_gene.get(gene, [])
    return (
        uniq_join(match.get("omim_phenotype", "") for match in matches),
        uniq_join(match.get("omim_inheritance", "") for match in matches),
    )


def join_orphanet(ensembl_gene_id_value, orphanet_by_ensg):
    return orphanet_by_ensg.get(ensembl_gene_id_value, {}).get("Orphanet", "0") or "0"


def join_imprinting(gene, imprinting_by_gene):
    match = imprinting_by_gene.get(gene, {})
    return match.get("Imprinting_status", ""), match.get("Imprinting_expressed_allele", "")


def join_pseudoautosomal(ensembl_gene_id_value, pseudo_by_ensg):
    return pseudo_by_ensg.get(ensembl_gene_id_value, {}).get("Pseudoautosomal", "")


def primary_constraint_values(transcript, constraint_by_transcript):
    keys = [
        "Gnomad_oe_lof_score",
        "Gnomad_oe_ci_lower",
        "Gnomad_oe_ci_upper",
        "Gnomad_oe_mis_score",
        "Gnomad_mis_z_score",
        "Gnomad_pLI_score",
        "Gnomad_pnull_score",
        "Gnomad_prec_score",
    ]
    if not present(transcript):
        return {k: "" for k in keys}
    transcript_no_version = transcript.split(".", 1)[0]
    match = constraint_by_transcript.get(transcript) or constraint_by_transcript.get(transcript_no_version) or {}
    return {k: match.get(k, "") for k in keys}


def constraint_all_summary(transcripts, constraint_by_transcript):
    parts = []
    for transcript in transcripts:
        if not present(transcript):
            continue
        values = primary_constraint_values(transcript, constraint_by_transcript)
        parts.append(
            f"{transcript}|lof={values['Gnomad_oe_lof_score']}|mis={values['Gnomad_oe_mis_score']}|pli={values['Gnomad_pLI_score']}"
        )
    return uniq_join(parts, sep=";")


def make_columns(mode, samples):
    sample_headers = [cre_sample_name(sample) for sample in samples]
    columns = ["Position", "UCSC_Link", "GNOMAD_Link", "Ref", "Alt"]
    columns.extend([f"Zygosity.{sample}" for sample in sample_headers])
    columns.append("Gene")
    columns.extend([f"Burden.{sample}" for sample in sample_headers])
    columns.extend(["gts", "Variation", "Info", "Refseq_change", "Depth", "Quality"])
    columns.extend([f"Alt_depths.{sample}" for sample in sample_headers])
    columns.extend([f"gt_quals.{sample}" for sample in sample_headers])
    if mode == "coding":
        columns.extend(
            [
                "Trio_coverage",
                "Ensembl_gene_id",
                "Gene_description",
                "omim_phenotype",
                "omim_inheritance",
                "Orphanet",
                "Clinvar",
                "HGMD_id",
                "HGMD_gene",
                "HGMD_tag",
                "HGMD_ref",
                "Gnomad_af_grpmax",
                "Gnomad_af",
                "Gnomad_ac",
                "Gnomad_hom",
                "Gnomad_male_ac",
                "Gnomad_fafmax_faf95_max",
                "Gnomad_filter",
                "Regeneron_exome_AF",
                "Regeneron_exome_AC",
                "Ensembl_transcript_id",
                "rsIDs",
                "AA_position",
                "Exon",
                "Protein_domains",
                "Gnomad_oe_lof_score",
                "Gnomad_oe_ci_lower",
                "Gnomad_oe_ci_upper",
                "Gnomad_oe_mis_score",
                "Gnomad_mis_z_score",
                "Gnomad_pLI_score",
                "Gnomad_pnull_score",
                "Gnomad_prec_score",
                "Sift_score",
                "Polyphen_score",
                "Cadd_score",
                "Vest4_score",
                "Revel_score",
                "Gerp_score",
                "AlphaMissense",
                "phylop100way",
                "SpliceAI_impact",
                "SpliceAI_score",
                "Imprinting_status",
                "Imprinting_expressed_allele",
                "Pseudoautosomal",
                "Old_multiallelic",
                "CSQ_biotype",
                "CSQ_impact",
                "Gene_all",
                "Ensembl_gene_id_all",
                "Ensembl_transcript_id_all",
                "Variation_all",
                "Info_all",
                "Refseq_change_all",
                "AA_position_all",
                "Exon_all",
                "Sift_score_all",
                "Polyphen_score_all",
                "Gene_description_all",
                "omim_phenotype_all",
                "omim_inheritance_all",
                "Orphanet_all",
                "constraint_all",
                "CSQ_biotype_all",
                "CSQ_impact_all",
            ]
        )
        return columns

    columns.extend(
        [
            "Trio_coverage",
            "Ensembl_gene_id",
            "Gene_description",
            "omim_phenotype",
            "omim_inheritance",
            "Orphanet",
            "Clinvar",
            "HGMD_id",
            "HGMD_gene",
            "HGMD_tag",
            "HGMD_ref",
            "Gnomad_af_grpmax",
            "Gnomad_af",
            "Gnomad_ac",
            "Gnomad_hom",
            "Gnomad_male_ac",
            "Gnomad_fafmax_faf95_max",
            "Gnomad_filter",
            "Regeneron_exome_AF",
            "Regeneron_exome_AC",
            "Ensembl_transcript_id",
            "rsIDs",
            "Gnomad_oe_lof_score",
            "Gnomad_oe_ci_lower",
            "Gnomad_oe_ci_upper",
            "Gnomad_oe_mis_score",
            "Gnomad_mis_z_score",
            "Gnomad_pLI_score",
            "Gnomad_pnull_score",
            "Gnomad_prec_score",
            "Cadd_score",
            "phylop100way",
            "SpliceAI_impact",
            "SpliceAI_score",
            "ncER_score",
            "ReMM_score",
            "LINSIGHT_score",
            "Noncoding_path_pred",
            "promoterAI_score",
            "Imprinting_status",
            "Imprinting_expressed_allele",
            "Pseudoautosomal",
            "Old_multiallelic",
            "TF_binding_sites",
            "GreenDB_variant_type",
            "GreenDB_closest_gene",
            "GreenDB_controlled_gene",
            "CTCF_binding_site",
            "DNaseI_hypersensitive_site",
            "ENH_cellline_tissue",
            "UCE_100bp",
            "UCE_200bp",
            "CSQ_biotype",
            "CSQ_impact",
            "Gene_all",
            "Ensembl_gene_id_all",
            "Ensembl_transcript_id_all",
            "Variation_all",
            "Info_all",
            "Refseq_change_all",
            "Gene_description_all",
            "omim_phenotype_all",
            "omim_inheritance_all",
            "Orphanet_all",
            "constraint_all",
            "CSQ_biotype_all",
            "CSQ_impact_all",
        ]
    )
    return columns


def set_row_value(row, field, value):
    if field in row:
        row[field] = value


def cre_sample_name(sample):
    return sample.replace("-", "_")


def consequence_display(csq, order_map):
    terms = consequence_terms(csq)
    if not terms:
        return ""
    return sorted(terms, key=lambda term: (order_map.get(normalize_consequence_term(term), 10**9), term))[0]


def split_hgvs_suffix(value):
    if not present(value):
        return "", ""
    text = str(value)
    if ":" in text:
        prefix, suffix = text.split(":", 1)
        return prefix, suffix
    return "", text


def preferred_refseq_accession(csq):
    for field in ("MANE_SELECT", "MANE_PLUS_CLINICAL"):
        value = mane_value(csq, field)
        if value.startswith(("NM_", "NR_", "XM_", "XR_")):
            return value
    feature = csq.get("Feature", "")
    if feature.startswith(("NM_", "NR_", "XM_", "XR_")):
        return feature
    return ""


def preferred_refseq_rank(csq):
    accession = preferred_refseq_accession(csq)
    if mane_select_value(csq):
        return (0, 0 if accession.startswith(("NM_", "NR_")) else 1)
    if mane_plus_clinical_value(csq):
        return (1, 0 if accession.startswith(("NM_", "NR_")) else 1)
    if accession.startswith(("NM_", "NR_")):
        return (2, 0)
    if accession.startswith(("XM_", "XR_")):
        return (3, 0)
    return (4, 0)


def refseq_change_for_csq(csq):
    accession = preferred_refseq_accession(csq)
    hgvsc_prefix, hgvsc_suffix = split_hgvs_suffix(csq.get("HGVSc", ""))
    _, hgvsp_suffix = split_hgvs_suffix(csq.get("HGVSp", ""))
    if not accession or not hgvsc_suffix:
        return ""
    if hgvsp_suffix:
        return f"{accession}:{hgvsc_suffix}:{hgvsp_suffix}"
    return f"{accession}:{hgvsc_suffix}"


def best_refseq_change_for_gene(csq_records, symbol, order_map, mode):
    candidates = []
    for csq in csq_records:
        if gene_symbol(csq) != symbol:
            continue
        refseq_change = refseq_change_for_csq(csq)
        if not present(refseq_change):
            continue
        candidates.append((preferred_refseq_rank(csq), csq_sort_key(csq, order_map, mode), refseq_change))
    if candidates:
        return sorted(candidates, key=lambda item: (item[0], item[1], item[2]))[0][2]
    return "NA"


def parse_args():
    parser = argparse.ArgumentParser(description="Unified slivar report builder for coding, wgs, and wgs-high-impact.")
    parser.add_argument("--mode", required=True, choices=["coding", "wgs", "wgs-high-impact"])
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--impact-order-file", required=True)
    parser.add_argument("--gene-descriptions", default="")
    parser.add_argument("--omim", default="")
    parser.add_argument("--orphanet", default="")
    parser.add_argument("--constraint", default="")
    parser.add_argument("--imprinting", default="")
    parser.add_argument("--pseudoautosomal", default="")
    return parser.parse_args()


def main():
    args = parse_args()
    order_map = load_impact_order(args.impact_order_file)
    gene_descriptions = index_first(load_table(args.gene_descriptions), "ensembl_gene_id")
    omim = index_many(load_table(args.omim), "gene_name")
    orphanet = index_first(load_table(args.orphanet), "Ensembl_gene_id")
    constraint = index_constraint(load_table(args.constraint))
    imprinting = index_first(load_table(args.imprinting), "Gene")
    pseudoautosomal = index_first(load_table(args.pseudoautosomal), "Ensembl_gene_id")

    rows = []
    with pysam.VariantFile(args.vcf) as vcf:
        samples = list(vcf.header.samples)
        sample_headers = {sample: cre_sample_name(sample) for sample in samples}
        columns = make_columns(args.mode, samples)
        csq_fields = parse_csq_header(vcf)
        if not csq_fields:
            raise SystemExit("VCF header does not contain a usable INFO/CSQ definition")
        grouped_records = group_records_by_variant(vcf)

        for key in sorted(grouped_records):
            records = grouped_records[key]
            record = records[0]
            csq_records = merged_csq_records(records, csq_fields)
            primary = choose_primary_csq(csq_records, order_map, args.mode)

            ref = record.ref
            alt = record.alts[0] if record.alts else ""
            position = f"{record.chrom}:{record.pos}"

            row = {column: "" for column in columns}
            row["Position"] = position
            row["UCSC_Link"] = hyperlink_ucsc(position)
            row["GNOMAD_Link"] = hyperlink_gnomad(record.chrom, record.pos, ref, alt)
            row["Ref"] = ref
            row["Alt"] = alt

            sample_gt_strings = []
            sample_depths = []
            sample_alt_depths = []
            sample_gqs = []
            for sample in samples:
                sample_data = record.samples[sample]
                sample_header = sample_headers[sample]
                sample_gt_strings.append(gt_string(sample_data, ref, record.alts or []))
                sample_depths.append(sample_depth(sample_data))
                sample_alt_depths.append(sample_alt_depth(sample_data))
                sample_gqs.append(sample_gq(sample_data))
                row[f"Zygosity.{sample_header}"] = zygosity(sample_data)
                row[f"Alt_depths.{sample_header}"] = sample_alt_depths[-1]
                row[f"gt_quals.{sample_header}"] = sample_gqs[-1]

            row["gts"] = ",".join(sample_gt_strings)
            row["Depth"] = as_text(first_info(records, "DP"))
            row["Quality"] = "" if record.qual is None else str(record.qual)
            row["Trio_coverage"] = "_".join(sample_depths)
            row["Clinvar"] = merged_clinvar_text(records) or "."
            row["HGMD_id"] = "NA"
            row["HGMD_tag"] = "NA"
            row["HGMD_ref"] = "NA"

            for field, source in [
                ("Gnomad_af_grpmax", "gnomad_af_grpmax"),
                ("Gnomad_af", "gnomad_af"),
                ("Gnomad_ac", "gnomad_ac"),
                ("Gnomad_hom", "gnomad_hom"),
                ("Gnomad_male_ac", "gnomad_male_ac"),
                ("Gnomad_fafmax_faf95_max", "gnomad_fafmax_faf95_max"),
                ("Gnomad_filter", "gnomad_filter"),
                ("Regeneron_exome_AF", "regeneron_exome_AF"),
                ("Regeneron_exome_AC", "regeneron_exome_AC"),
                ("rsIDs", "rs_ids"),
                ("Cadd_score", "CADD_phred"),
                ("Vest4_score", "Vest4_score"),
                ("Revel_score", "REVEL_score"),
                ("Gerp_score", "Gerp_score"),
                ("AlphaMissense", "AlphaMissense"),
                ("phylop100way", "phyloP100way"),
                ("ncER_score", "ncER"),
                ("ReMM_score", "ReMM"),
                ("LINSIGHT_score", "LinSight_Score"),
                ("TF_binding_sites", "tf_binding_sites"),
                ("GreenDB_variant_type", "GreenDB_variant_type"),
                ("GreenDB_closest_gene", "GreenDB_closest_gene"),
                ("GreenDB_controlled_gene", "GreenDB_controlled_gene"),
                ("CTCF_binding_site", "ctcf_binding_site"),
                ("DNaseI_hypersensitive_site", "dnasei_hypersensitive_site"),
                ("ENH_cellline_tissue", "enh_cellline_tissue"),
                ("UCE_100bp", "uce_100bp"),
                ("UCE_200bp", "uce_200bp"),
                ("Old_multiallelic", "OLD_MULTIALLELIC"),
            ]:
                set_row_value(row, field, as_text(first_info(records, source)))

            if "Gnomad_filter" in row:
                row["Gnomad_filter"] = row["Gnomad_filter"] or "None"
            if "Cadd_score" in row:
                row["Cadd_score"] = row["Cadd_score"] or "None"
            for field in ["Vest4_score", "Revel_score", "Gerp_score", "AlphaMissense", "ncER_score", "ReMM_score", "LINSIGHT_score"]:
                if field in row:
                    row[field] = row.get(field, "") or "None"

            spliceai_impact, spliceai_score = parse_spliceai(first_info(records, "spliceai_score"))
            set_row_value(row, "SpliceAI_impact", spliceai_impact)
            set_row_value(row, "SpliceAI_score", spliceai_score)
            set_row_value(row, "promoterAI_score", parse_promoterai(first_info(records, "promoterAI")))
            if "Noncoding_path_pred" in row:
                row["Noncoding_path_pred"] = noncoding_pred_fraction(
                    row.get("Cadd_score", ""), row.get("ncER_score", ""), row.get("ReMM_score", ""), row.get("LINSIGHT_score", "")
                )

            all_gene_symbols = [gene_symbol(csq) for csq in csq_records if present(gene_symbol(csq))]
            all_ensg = [best_ensembl_gene_id_for_symbol(csq_records, symbol) for symbol in all_gene_symbols]
            all_ensg = [g for g in all_ensg if present(g)]
            all_ensembl_transcripts = [csq.get("Feature", "") for csq in csq_records if csq.get("Feature", "").startswith("ENST")]

            row["Gene_all"] = uniq_join(all_gene_symbols)
            row["Ensembl_gene_id_all"] = uniq_join(all_ensg)
            row["Ensembl_transcript_id_all"] = uniq_join(all_ensembl_transcripts)
            row["Variation_all"] = uniq_join(csq.get("Consequence", "") for csq in csq_records)
            row["Info_all"] = build_info_all(csq_records)
            row["Refseq_change_all"] = refseq_change_all(csq_records)
            set_row_value(row, "AA_position_all", uniq_join(clean_slash(csq.get("Protein_position", "")) for csq in csq_records))
            set_row_value(row, "Exon_all", uniq_join(clean_slash(csq.get("EXON", "")) for csq in csq_records))
            set_row_value(row, "Sift_score_all", uniq_join(parse_vep_score(csq.get("SIFT", "")) for csq in csq_records) or "None")
            set_row_value(row, "Polyphen_score_all", uniq_join(parse_vep_score(csq.get("PolyPhen", "")) for csq in csq_records) or "None")
            row["CSQ_biotype_all"] = uniq_join(csq.get("BIOTYPE", "") for csq in csq_records)
            row["CSQ_impact_all"] = uniq_join(csq.get("IMPACT", "") for csq in csq_records)

            gene_desc_all = []
            omim_pheno_all = []
            omim_inh_all = []
            orphanet_all = []
            for symbol in row["Gene_all"].split(","):
                if not present(symbol):
                    continue
                ensg = best_ensembl_gene_id_for_symbol(csq_records, symbol)
                if present(ensg):
                    gene_desc_all.append(join_gene_description(ensg, gene_descriptions))
                    orphanet_all.append(join_orphanet(ensg, orphanet))
                phen, inh = join_omim(symbol, omim)
                omim_pheno_all.append(phen)
                omim_inh_all.append(inh)
            row["Gene_description_all"] = uniq_join(gene_desc_all)
            row["omim_phenotype_all"] = uniq_join(omim_pheno_all)
            row["omim_inheritance_all"] = uniq_join(omim_inh_all)
            row["Orphanet_all"] = uniq_join(orphanet_all)
            row["constraint_all"] = constraint_all_summary(row["Ensembl_transcript_id_all"].split(","), constraint)

            if primary is not None:
                primary_gene = gene_symbol(primary)
                primary_ensg = best_ensembl_gene_id_for_symbol(csq_records, primary_gene)
                primary_feature = primary.get("Feature", "")
                primary_tx = primary_feature if primary_feature.startswith("ENST") else ""
                row["Gene"] = primary_gene
                row["Ensembl_gene_id"] = primary_ensg
                row["Variation"] = consequence_display(primary, order_map)
                row["Info"] = build_info_item(primary) or "NA"
                row["Refseq_change"] = best_refseq_change_for_gene(csq_records, primary_gene, order_map, args.mode)
                row["Ensembl_transcript_id"] = primary_tx
                set_row_value(row, "AA_position", clean_slash(primary.get("Protein_position", "")))
                set_row_value(row, "Exon", clean_slash(primary.get("EXON", "")))
                set_row_value(row, "Protein_domains", primary.get("DOMAINS", ""))
                set_row_value(row, "Sift_score", parse_vep_score(primary.get("SIFT", "")) or "None")
                set_row_value(row, "Polyphen_score", parse_vep_score(primary.get("PolyPhen", "")) or "None")
                row["CSQ_biotype"] = primary.get("BIOTYPE", "")
                row["CSQ_impact"] = primary.get("IMPACT", "")
                row["Gene_description"] = join_gene_description(primary_ensg, gene_descriptions)
                row["omim_phenotype"], row["omim_inheritance"] = join_omim(primary_gene, omim)
                row["Orphanet"] = join_orphanet(primary_ensg, orphanet)
                row["Imprinting_status"], row["Imprinting_expressed_allele"] = join_imprinting(primary_gene, imprinting)
                row["Pseudoautosomal"] = join_pseudoautosomal(primary_ensg, pseudoautosomal)
                row["HGMD_gene"] = primary_gene or "NA"
                for field, value in primary_constraint_values(primary_tx, constraint).items():
                    row[field] = value
            else:
                row["Info"] = "NA"
                row["Refseq_change"] = "NA"
                row["Orphanet"] = "0"
                row["HGMD_gene"] = "NA"

            rows.append(row)

    burden_seen = {sample: set() for sample in samples}
    burden_by_sample = {sample: Counter() for sample in samples}
    for row in rows:
        gene = row["Gene"]
        if not present(gene):
            continue
        variant_gene_key = (row["Position"], row["Ref"], row["Alt"], gene)
        for sample in samples:
            if row[f"Zygosity.{sample_headers[sample]}"] in {"Het", "Hom"} and variant_gene_key not in burden_seen[sample]:
                burden_seen[sample].add(variant_gene_key)
                burden_by_sample[sample][gene] += 1
    for row in rows:
        for sample in samples:
            row[f"Burden.{sample_headers[sample]}"] = str(burden_by_sample[sample].get(row["Gene"], 0))

    rows.sort(key=lambda row: row["Position"])
    with open(args.out_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
