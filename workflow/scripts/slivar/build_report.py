#!/usr/bin/env python3

import argparse
import csv
import os
from collections import Counter, defaultdict

import pysam

from shared import (
    MISSING,
    find_ensembl_gene_id,
    format_genotype,
    get_alt_depth,
    get_csq_ensembl_gene_id,
    get_csq_fields,
    get_gene_symbol,
    get_info_value,
    get_sample_depth,
    get_sample_gq,
    has_value,
    load_consequence_order,
    parse_csq_annotations,
    primary_csq_sort_key,
    select_consequence,
    select_ensembl_gene_id,
    select_primary_csq,
    value_as_text,
)


# Transcript constraint columns used by both primary and *_all report fields.
CONSTRAINT_COLUMNS = [
    "Gnomad_oe_lof_score",
    "Gnomad_oe_ci_lower",
    "Gnomad_oe_ci_upper",
    "Gnomad_oe_mis_score",
    "Gnomad_mis_z_score",
    "Gnomad_pLI_score",
    "Gnomad_pnull_score",
    "Gnomad_prec_score",
]


# Fields where empty report values are written as "." to match CRE output.
DOT_MISSING_FIELDS = {
    "AA_position", "AlphaMissense", "Cadd_score", "Clinvar", "ENH_cellline_tissue",
    # Ensembl_gene_id is deliberately absent: annotate_compound_hets.py left-merges the
    # compound-het table onto the report on this column, and compound_hets.py groups
    # gene-less variants under the literal ".". Writing "." here would make every
    # gene-less row match that pseudo-gene and inherit its CH_status. Leave it empty --
    # Gemini does the same, and annotate_compound_hets fills it back to "." on output.
    "Ensembl_gene_id_all", "Ensembl_transcript_id", "Exon", "Gene", "Gene_all",
    "Gene_description_all", "Gerp_score", "Gnomad_af", "Gnomad_af_grpmax",
    "Gnomad_fafmax_faf95_max", "Gnomad_filter", "Gnomad_hom", "Gnomad_male_ac",
    "Gnomad_mis_z_score", "Gnomad_oe_ci_lower", "Gnomad_oe_ci_upper",
    "Gnomad_oe_lof_score", "Gnomad_oe_mis_score", "Gnomad_pLI_score",
    "Gnomad_pnull_score", "Gnomad_prec_score", "GSO_AC", "GSO_AF", "GSO_hemi",
    "GSO_nhomalt", "GreenDB_closest_gene", "GreenDB_controlled_gene",
    "GreenDB_variant_type", "HGMD_gene", "HGMD_id", "HGMD_ref", "HGMD_tag",
    "Imprinting_expressed_allele", "Imprinting_status", "LINSIGHT_score",
    "Old_multiallelic", "Orphanet_all", "Polyphen_score", "Polyphen_score_all",
    "Protein_domains", "Pseudoautosomal", "Quality", "ReMM_score", "Refseq_change",
    "Refseq_change_all", "Revel_score", "Sift_score", "Sift_score_all",
    "TF_binding_sites", "Vest4_score", "ncER_score", "omim_inheritance_all",
    "omim_phenotype_all", "promoterAI_score", "phylop100way", "rsIDs",
    "CoLoRSdb_AF", "CoLoRSdb_AC", "CoLoRSdb_AC_Hemi", "CoLoRSdb_nhomalt",
    "TG_LRWGS_samples", "Dark_genes", "PS",
}


# Numeric frequency/count fields where missing values are reported as 0.
ZERO_MISSING_FIELDS = {
    "GSO_AC", "GSO_AF", "GSO_hemi", "GSO_nhomalt",
    "CoLoRSdb_AF", "CoLoRSdb_AC", "CoLoRSdb_AC_Hemi", "CoLoRSdb_nhomalt",
    "Gnomad_ac", "Gnomad_af", "Gnomad_af_grpmax", "Gnomad_fafmax_faf95_max",
    "Gnomad_hom", "Regeneron_exome_AF", "Regeneron_exome_AC",
    "thousandG_AF", "thousandG_AC", "thousandG_nhomalt",
    "TG_LRWGS_AC", "TG_LRWGS_hom",
}


# INFO flags become explicit 1/0 report columns.
FLAG_BINARY_FIELDS = {"CTCF_binding_site", "DNaseI_hypersensitive_site", "UCE_100bp", "UCE_200bp"}


# GreenDB INFO values must retain empty entries so overlapping-region fields stay aligned.
ALIGNED_LIST_FIELDS = {"GreenDB_variant_type", "GreenDB_closest_gene", "GreenDB_controlled_gene"}


# Sample alt-depth columns also use 0 for missing values.
ZERO_MISSING_PREFIXES = ("Alt_depths.",)


# Helper for turning raw annotation/reference inputs into report strings.
def join_unique_values(values, sep=","):
    # seen handles duplicate checks; output keeps the first-seen report order.
    seen = set()
    output = []
    for value in values:
        if value is None:
            continue
        # Some annotations arrive as tuple/list values; wrap scalars so both paths iterate the same way.
        items = value if isinstance(value, (tuple, list)) else [value]
        for item in items:
            text = str(item)
            # Drop standard missing tokens and values already added for this field.
            if text in MISSING or text in seen:
                continue
            seen.add(text)
            output.append(text)
    return sep.join(output)


def replace_slashes(value):
    if value is None:
        return ""
    return str(value).replace("/", "_")


# Read a CSV/TSV reference table as one dict per row.
# If delimiter is not provided, infer comma for .csv and tab otherwise.
def load_reference_table(path, delimiter=None):
    if not path:
        return []
    if delimiter is None:
        delimiter = "," if path.endswith(".csv") else "\t"
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


# Index tables by key, keeping the first row for duplicate keys.
def index_first_by(rows, key_column):
    index = {}
    for row in rows:
        key = row.get(key_column, "")
        if key and key not in index:
            index[key] = row
    return index


# Index tables by key, keeping every row for joins such as OMIM.
def index_many_by(rows, key_column):
    index = defaultdict(list)
    for row in rows:
        key = row.get(key_column, "")
        if key:
            index[key].append(row)
    return index

# Report columns replace hyphens in sample IDs with underscores.
def report_sample_name(sample):
    return sample.replace("-", "_")


# Convert one sample's FORMAT fields into the zygosity label.
def classify_zygosity(sample_data, ref, alts, chrom=""):
    # format_genotype changes GT allele indexes into allele text, e.g. 0/1 -> A/G.
    genotype = format_genotype(sample_data, ref, alts).replace("|", "/").replace("./.", "Missing")
    if "Missing" in genotype:
        return "Missing"
    # A called genotype with no ALT read support is reported as reference/no-call.
    if get_alt_depth(sample_data) == 0:
        return "-"
    alleles = genotype.split("/")
    if len(alleles) == 2:
        # Diploid examples for REF=A, ALT=G: A/A -> "-", G/G -> "Hom", A/G -> "Het".
        if alleles[0] == alleles[1]:
            return "-" if alleles[0] == ref else "Hom"
        return "Het"
    chrom_text = str(chrom)
    if "X" in chrom_text or "Y" in chrom_text:
        # Haploid X/Y calls can arrive as one allele: "." missing, REF "-", ALT "Hom".
        if genotype == ".":
            return "Missing"
        return "-" if genotype == ref else "Hom"
    return genotype


# Annotation-specific formatting used while constructing each report row.
def parse_vep_score(value):
    if not has_value(value):
        return ""
    text = str(value)
    if "(" in text and ")" in text:
        return text.rsplit("(", 1)[1].split(")", 1)[0]
    return text

# The GreenDB columns hold one entry per overlapping region and are read across
# by position, so an entry with no value has to keep its slot rather than vanish.
def as_aligned_text(value):
    if not has_value(value):
        return ""
    items = value if isinstance(value, tuple) else (value,)
    return ",".join("" if item is None or str(item) in MISSING else str(item) for item in items)


def format_phase_sets(value):
    if isinstance(value, tuple):
        return "|".join(str(item) if has_value(item) else "." for item in value)
    return value_as_text(value).replace(",", "|")


# Apply missing-value conventions after all row fields are populated.
def normalize_report_value(field, value):
    text = "" if value is None else str(value).strip()
    if (field in ZERO_MISSING_FIELDS or field.startswith(ZERO_MISSING_PREFIXES)) and text in {"", ".", "-1", "NA", "None"}:
        return "0"
    if field in DOT_MISSING_FIELDS and text in {"", "NA", "None"}:
        return "."
    if field in FLAG_BINARY_FIELDS:
        return "1" if text == "1" else "0"
    return value

# Build one transcript summary for Info/Info_all:
# Gene:exon<EXON>:HGVSc:HGVSp.
def build_info_item(csq):
    gene = get_gene_symbol(csq)
    # Missing exon text is reported as exonNA rather than an empty exon label.
    exon = csq.get("EXON", "") or "NA"
    if not has_value(exon):
        exon = "NA"
    # VEP percent-encodes "=" as "%3D"; reports use the readable character.
    return f"{gene}:exon{exon}:{csq.get('HGVSc', '')}:{csq.get('HGVSp', '')}".replace("%3D", "=")


# Keep only the c./p. suffix so it can be paired with the chosen RefSeq accession.
def split_hgvs_suffix(value):
    if not has_value(value):
        return "", ""
    text = str(value)
    if ":" in text:
        prefix, suffix = text.split(":", 1)
        return prefix, suffix
    return "", text


# Return the RefSeq accession available on this CSQ.
def preferred_refseq_accession(csq):
    # Prefer MANE RefSeq accessions, then fall back to a RefSeq-like Feature value.
    for field in ("MANE_SELECT", "MANE_PLUS_CLINICAL"):
        value = csq.get(field, "")
        if not has_value(value):
            value = ""
        if value.startswith(("NM_", "NR_", "XM_", "XR_")):
            return value
    feature = csq.get("Feature", "")
    if feature.startswith(("NM_", "NR_", "XM_", "XR_")):
        return feature
    return ""

# Lower tuple values win: MANE, curated NM/NR accessions, then predicted XM/XR.
def preferred_refseq_rank(csq):
    accession = preferred_refseq_accession(csq)
    if has_value(csq.get("MANE_SELECT", "")):
        source_rank = 0
    elif has_value(csq.get("MANE_PLUS_CLINICAL", "")):
        source_rank = 1
    elif accession.startswith(("NM_", "NR_")):
        return (2, 0)
    elif accession.startswith(("XM_", "XR_")):
        return (3, 0)
    else:
        return (4, 0)

    # For MANE ties, prefer curated NM_/NR_ accessions over predicted XM_/XR_.
    accession_rank = 0 if accession.startswith(("NM_", "NR_")) else 1
    return (source_rank, accession_rank)


def refseq_change_for_csq(csq):
    # Build accession:c.change[:p.change]; skip CSQs without a RefSeq accession or HGVSc.
    accession = preferred_refseq_accession(csq)
    _, hgvsc_suffix = split_hgvs_suffix(csq.get("HGVSc", ""))
    _, hgvsp_suffix = split_hgvs_suffix(csq.get("HGVSp", ""))
    if not accession or not hgvsc_suffix:
        return ""
    if hgvsp_suffix:
        return f"{accession}:{hgvsc_suffix}:{hgvsp_suffix}".replace("%3D", "=")
    return f"{accession}:{hgvsc_suffix}".replace("%3D", "=")

# Build the *_all field from every CSQ that can form a RefSeq HGVS string.
def refseq_change_all(csq_records):
    refseq_changes = []
    for csq in csq_records:
        refseq_change = refseq_change_for_csq(csq)
        # Empty CSQs are skipped; join_unique_values removes repeats but keeps first-seen order.
        if has_value(refseq_change):
            refseq_changes.append(refseq_change)
    return join_unique_values(refseq_changes) or "NA"

# Refseq_change is RefSeq-specific: keep it on the reported Gene, but allow a
# different CSQ than Info so an Ensembl-primary call can still report the best RefSeq change.
def best_refseq_change_for_gene(csq_records, symbol, order_map):
    candidates = []
    for csq in csq_records:
        if get_gene_symbol(csq) != symbol:
            continue
        refseq_change = refseq_change_for_csq(csq)
        if not has_value(refseq_change):
            continue
        # Sort by RefSeq quality, then report-CSQ preference, then text for stable ties.
        candidates.append((preferred_refseq_rank(csq), primary_csq_sort_key(csq, order_map), refseq_change))
    if candidates:
        return sorted(candidates, key=lambda item: (item[0], item[1], item[2]))[0][2]
    return "NA"


# Build the report SpliceAI_impact and SpliceAI_score fields.
def parse_spliceai(value):
    if not has_value(value):
        # No SpliceAI record for this variant. A real max score of 0 is reported as 0 below.
        return "NA|NA|NA", "."
    # Roughly follows the SpliceAI parsing used by postfilter.py and the former CRE/R implementation.
    # Values look like Allele|Gene|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|...
    # Repeated INFO values might be a tuple. Each value can also still hold
    # comma-separated annotations, so flatten both forms first.
    annotations = []
    if isinstance(value, tuple):
        for item in value:
            annotations.extend(str(item).split(","))
    else:
        annotations = str(value).split(",")
    # These report fields stay tied together: the highest DS score determines
    # which gene, splice impact label, and DP position are reported.
    best_gene = "NA"
    best_impact = "NA"
    best_score = 0.0
    best_pos = "NA"
    fields = [("acceptor_gain", 2, 6), ("acceptor_loss", 3, 7), ("donor_gain", 4, 8), ("donor_loss", 5, 9)]
    # Track the strongest DS_* value along with its gene, impact type, and DP position.
    for annotation in annotations:
        parts = annotation.split("|")
        # Short or malformed annotations cannot contain all four DS/DP pairs.
        if len(parts) < 10:
            continue
        scores = []
        gene = parts[1]
        for impact_name, score_idx, pos_idx in fields:
            try:
                score = float(parts[score_idx])
            except Exception:
                continue
            scores.append((impact_name, score, parts[pos_idx]))
        if not scores:
            continue
        # Pick the strongest DS value within this annotation, keeping its impact and position.
        max_score = max(score for _, score, _ in scores)
        impact = "NA"
        pos = "NA"
        for impact_name, score, pos_value in scores:
            if score == max_score:
                impact = impact_name
                pos = pos_value
        # Only replace the current report values when this annotation has the
        # strongest score seen so far.
        if best_score < max_score:
            best_score = max_score
            best_gene = gene
            best_impact = impact
            best_pos = pos
    if best_score == 0:
        return "NA|NA|NA", "0"
    return f"{best_gene}|{best_impact}|{best_pos}", str(best_score)


def max_numeric_csv(value):
    # Some INFO fields contain comma-separated numeric predictions; report the largest.
    if not has_value(value):
        return ""
    best_score = None
    best_text = ""
    items = value if isinstance(value, tuple) else [value]
    for item in items:
        for text in str(item).split(","):
            text = text.strip()
            try:
                score = float(text)
            except Exception:
                continue
            if best_score is None or score > best_score:
                best_score = score
                best_text = text
    if best_score is None:
        return ""
    return best_text


# Build the report promoterAI_score field while preserving the score text.
def parse_promoterai(value):
    if not has_value(value):
        return "."
    # Roughly follows the promoterAI parsing used by postfilter.py and the former CRE/R implementation.
    # Values are numeric score strings, for example "0.22" or "-0.14,0.03".
    # Repeated INFO values might be a tuple; otherwise split comma-separated scores.
    raw_values = value if isinstance(value, tuple) else str(value).split(",")
    kept = []
    for raw in raw_values:
        text = str(raw).strip()
        # Missing values and "No" promoterAI hits are omitted; no kept scores reports as ".".
        if text in MISSING or text in {"No", "no"}:
            continue
        kept.append(text)
    return ",".join(kept) if kept else "."


# Match the former CRE/R noncoding_pred function: report passing scores as "pass//available".
def noncoding_pred_fraction(cadd, ncer, remm, linsight):
    # Missing/non-numeric scores are ignored and do not count in the denominator.
    def float_or_none(value):
        if not has_value(value):
            return None
        try:
            return float(str(value))
        except Exception:
            return None

    preds = []
    # Thresholds are the former CRE/R cutoffs for CADD, ncER, ReMM, and LINSIGHT.
    values = [
        (float_or_none(cadd), 10),
        (float_or_none(ncer), 95.95),
        (float_or_none(remm), 0.9585),
        (float_or_none(linsight), 0.9828),
    ]
    for value, threshold in values:
        if value is not None:
            # True means this usable score predicts a noncoding impact.
            preds.append(value > threshold)
    return "0//0" if not preds else f"{sum(preds)}//{len(preds)}"


def hyperlink_ucsc(position):
    return '=HYPERLINK("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hgt.out3=10x&position=' + position + '","UCSC_link")'


def hyperlink_gnomad(chrom, pos, ref, alt):
    variant_id = f"{chrom}-{pos}-{ref}-{alt}"
    return '=HYPERLINK("http://gnomad.broadinstitute.org/variant/' + variant_id + '?dataset=gnomad_r4","GNOMAD_link")'

# Combine the ClinVar INFO fields while preserving first-seen unique text.
def clinvar_text(record):
    values = []
    for field in ("clinvar_pathogenic", "clinvar_sig", "clinvar_sig_conf"):
        value = get_info_value(record, field)
        if not has_value(value):
            continue
        if isinstance(value, tuple):
            for item in value:
                if has_value(item):
                    values.append(str(item))
        else:
            values.append(str(value))
    return join_unique_values(values, sep=";")


# OMIM is keyed by gene symbol and can have multiple phenotype rows per gene.
def join_omim(gene, omim_by_gene):
    matches = omim_by_gene.get(gene, [])
    phenotypes = []
    inheritances = []
    for match in matches:
        # Keep the two OMIM report columns separate while collecting every matching row.
        phenotypes.append(match.get("omim_phenotype", ""))
        inheritances.append(match.get("omim_inheritance", ""))
    # Drop blanks/repeats and keep first-seen order for each OMIM column.
    return (
        join_unique_values(phenotypes),
        join_unique_values(inheritances),
    )

# The Orphanet table writes 0 for a gene with no entry; report that as missing.
def lookup_orphanet(ensembl_gene_id_value, orphanet_by_ensg):
    value = orphanet_by_ensg.get(ensembl_gene_id_value, {}).get("Orphanet", "")
    if str(value).strip() in {"0", ""}:
        return ""
    return value

# Build the ordered gene-symbol list used for gene-level report lookups.
def genes_primary_first(primary_gene, gene_all_text):
    genes = [primary_gene] if has_value(primary_gene) else []
    for gene in gene_all_text.split(","):
        if has_value(gene) and gene not in genes:
            genes.append(gene)
    return genes


# Imprinting is symbol-keyed; check all report genes because the imprinted gene
# at a locus is not always the same gene chosen as primary.
def lookup_imprinting(genes, imprinting_by_gene):
    for gene in genes:
        match = imprinting_by_gene.get(gene)
        if match:
            return match.get("Imprinting_status", ""), match.get("Imprinting_expressed_allele", "")
    return "", ""


# HGMD_gene is a gene-level flag, so return the first report gene found in HGMD.
def first_hgmd_gene(genes, hgmd_genes):
    for gene in genes:
        if gene in hgmd_genes:
            return gene
    return "NA"


# Pseudoautosomal is keyed by Ensembl gene ID rather than display gene symbol.
def lookup_pseudoautosomal(ensembl_gene_ids, pseudo_by_ensg):
    for gene_id in ensembl_gene_ids:
        value = pseudo_by_ensg.get(gene_id, {}).get("Pseudoautosomal", "")
        if has_value(value):
            return value
    return ""


def load_hgmd(path):
    # HGMD is keyed by chrom:pos-ref-alt for exact variant-level matches.
    by_variant = {}
    # Also keep every HGMD gene symbol so HGMD_gene can be filled from the
    # selected/overlapping report genes even when the exact variant match is absent.
    genes = set()
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        for row in reader:
            # Required fields are chrom, pos, id, ref, alt, gene, and tag.
            if len(row) < 7:
                continue
            row = [value.strip() for value in row]
            if row[0].lower() == "chrom":
                continue
            # Citation fields are optional; pad them so the unpack below is stable.
            row.extend([""] * (13 - len(row)))
            chrom, pos, hgmd_id, ref, alt, hgmd_gene, hgmd_tag, author, allname, vol, page, year, pmid = row[:13]

            if has_value(hgmd_gene):
                genes.add(hgmd_gene)
            # Skip rows without the fields needed to build an exact variant key.
            if not all(has_value(value) for value in (chrom, pos, ref, alt)):
                continue

            variant_key_text = f"{chrom}:{pos}-{ref}-{alt}"
            # Start the per-variant HGMD bucket the first time this key is seen.
            match = by_variant.setdefault(
                variant_key_text,
                {"HGMD_id": [], "HGMD_tag": [], "HGMD_ref": []},
            )
            # Add this row's HGMD values; repeated rows are collapsed below.
            match["HGMD_id"].append(hgmd_id)
            match["HGMD_tag"].append(hgmd_tag)

            reference_parts = [author, allname, vol, page, year, "PMID:", pmid]
            # Add a citation string only when at least one reference field is present.
            if any(has_value(value) for value in [author, allname, vol, page, year, pmid]):
                match["HGMD_ref"].append(" ".join(reference_parts).strip())

    # Collapse repeated HGMD matches into the semicolon-delimited report fields.
    final_by_variant = {}
    for key, fields in by_variant.items():
        final_by_variant[key] = {}
        for field, values in fields.items():
            joined_values = join_unique_values(values, sep=";")
            final_by_variant[key][field] = joined_values or "NA"

    return final_by_variant, genes

# Return exactly the constraint columns used by the report for one transcript.
def primary_constraint_values(transcript, constraint_by_transcript):
    if not has_value(transcript):
        return {column: "" for column in CONSTRAINT_COLUMNS}
    transcript_no_version = transcript.split(".", 1)[0]
    match = constraint_by_transcript.get(transcript) or constraint_by_transcript.get(transcript_no_version) or {}
    return {column: match.get(column, "") for column in CONSTRAINT_COLUMNS}

# Compact transcript-level constraint values into the CRE *_all summary format.
def constraint_all_summary(transcripts, constraint_by_transcript):
    parts = []
    for transcript in transcripts:
        if not has_value(transcript):
            continue
        values = primary_constraint_values(transcript, constraint_by_transcript)
        parts.append(
            f"{transcript}|lof={values['Gnomad_oe_lof_score']}|mis={values['Gnomad_oe_mis_score']}|pli={values['Gnomad_pLI_score']}"
        )
    return join_unique_values(parts, sep=";")

# Build the exact report header sequence for the selected Slivar mode.
def make_columns(mode, samples, include_denovo=False, include_denovo_quality=False, profile="dragen"):
    sample_headers = [report_sample_name(sample) for sample in samples]
    # Each extend appends the next CSV header group; changing this order changes the report.
    columns = ["Position", "UCSC_Link", "GNOMAD_Link", "Ref", "Alt"]
    columns.extend([f"Zygosity.{sample}" for sample in sample_headers])
    columns.extend(["Gene", "Gene_all"])
    columns.extend([f"Burden.{sample}" for sample in sample_headers])
    columns.append("gts")
    if profile == "pacbio":
        columns.append("PS")
    columns.extend(["Variation", "Info", "Refseq_change", "Depth", "Quality"])
    columns.extend([f"Alt_depths.{sample}" for sample in sample_headers])
    columns.extend([f"gt_quals.{sample}" for sample in sample_headers])
    if include_denovo:
        columns.extend([f"denovo.{sample}" for sample in sample_headers])
    if include_denovo_quality:
        columns.extend([f"denovo_quality.{sample}" for sample in sample_headers])

    # Shared annotation columns used by both coding and WGS-style reports.
    columns.extend(["Trio_coverage", "Ensembl_gene_id", "Ensembl_gene_id_all"])
    columns.extend(["Gene_description_all", "omim_phenotype_all", "omim_inheritance_all", "Orphanet_all"])
    columns.extend(["Clinvar", "HGMD_id", "HGMD_gene", "HGMD_tag", "HGMD_ref"])
    columns.extend(["Gnomad_af_grpmax", "Gnomad_af", "Gnomad_ac", "Gnomad_hom", "Gnomad_male_ac"])
    columns.extend(["Gnomad_fafmax_faf95_max", "Gnomad_filter", "Regeneron_exome_AF", "Regeneron_exome_AC"])
    if profile == "pacbio":
        columns.extend(["CoLoRSdb_AF", "CoLoRSdb_AC", "CoLoRSdb_AC_Hemi", "CoLoRSdb_nhomalt"])
        columns.extend(["TG_LRWGS_AC", "TG_LRWGS_samples", "TG_LRWGS_hom"])
    else:
        columns.extend(["thousandG_AF", "thousandG_AC", "thousandG_nhomalt"])
        columns.extend(["GSO_AF", "GSO_AC", "GSO_nhomalt", "GSO_hemi"])
    columns.extend(["Ensembl_transcript_id", "rsIDs"])

    if mode == "coding":
        # Coding reports use coding-specific consequence and protein columns.
        columns.extend(["AA_position", "Exon", "Protein_domains"])
        columns.extend(CONSTRAINT_COLUMNS)
        columns.extend(["Sift_score", "Polyphen_score", "Cadd_score", "Vest4_score", "Revel_score", "Gerp_score", "AlphaMissense"])
        columns.extend(["phylop100way", "SpliceAI_impact", "SpliceAI_score"])
        columns.extend(["Imprinting_status", "Imprinting_expressed_allele", "Pseudoautosomal", "Old_multiallelic"])
        if profile == "pacbio":
            columns.append("Dark_genes")
        columns.extend(["CSQ_biotype", "CSQ_impact", "Sift_score_all", "Polyphen_score_all", "CSQ_biotype_all"])
    else:
        # All other Slivar modes use the WGS-style noncoding/report column set.
        columns.extend(CONSTRAINT_COLUMNS)
        columns.extend(["Cadd_score", "phylop100way", "SpliceAI_impact", "SpliceAI_score"])
        columns.extend(["ncER_score", "ReMM_score", "LINSIGHT_score", "Noncoding_path_pred", "promoterAI_score"])
        columns.extend(["Imprinting_status", "Imprinting_expressed_allele", "Pseudoautosomal", "Old_multiallelic"])
        if profile == "pacbio":
            columns.append("Dark_genes")
        columns.extend(["UCE_100bp", "UCE_200bp", "DNaseI_hypersensitive_site", "CTCF_binding_site"])
        columns.extend(["ENH_cellline_tissue", "TF_binding_sites"])
        if profile != "pacbio":
            columns.extend(["GreenDB_variant_type", "GreenDB_closest_gene", "GreenDB_controlled_gene"])
        columns.extend(["CSQ_biotype", "CSQ_impact", "CSQ_biotype_all"])

    # Put these aggregated annotation columns at the end of both report types.
    columns.extend(["Info_all", "CSQ_impact_all", "Ensembl_transcript_id_all", "Variation_all", "Refseq_change_all", "Constraint_all"])
    return columns


# Shared row-building code touches mode-specific columns; skip fields absent from this report type.
def set_if_column_present(row, column, value):
    if column in row:
        row[column] = value

# Drop optional de novo columns only when the source FORMAT fields were empty.
def drop_empty_optional_columns(columns, rows, prefixes):
    drop_columns = {
        column
        for column in columns
        if column.startswith(prefixes) and not any(has_value(row.get(column, "")) for row in rows)
    }
    if not drop_columns:
        return columns
    for row in rows:
        for column in drop_columns:
            row.pop(column, None)
    return [column for column in columns if column not in drop_columns]


def parse_args():
    parser = argparse.ArgumentParser(description="Unified slivar report builder.")
    parser.add_argument(
        "--mode",
        required=True,
        choices=["coding", "wgs-high-impact", "denovo", "panel", "panel-flank"],
    )
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--impact-order-file", required=True)
    parser.add_argument("--slivar-data-dir", required=True)
    parser.add_argument("--hgmd", default="")
    parser.add_argument("--profile", choices=["dragen", "pacbio"], default="dragen")
    return parser.parse_args()


def main():
    args = parse_args()
    order_map = load_consequence_order(args.impact_order_file)
    # Load the local reference tables used to fill report columns.
    gene_descriptions = index_first_by(load_reference_table(os.path.join(args.slivar_data_dir, "ensembl_w_description.txt")), "ensembl_gene_id")
    omim = index_many_by(load_reference_table(os.path.join(args.slivar_data_dir, "OMIM_hgnc_join_omim_phenos_2026-06-02.tsv")), "gene_name")
    orphanet = index_first_by(load_reference_table(os.path.join(args.slivar_data_dir, "orphanet.txt")), "Ensembl_gene_id")
    constraint = index_first_by(load_reference_table(os.path.join(args.slivar_data_dir, "gnomad_scores_transcript_level_v4.1.1.csv")), "Ensembl_transcript_id")
    imprinting = index_first_by(load_reference_table(os.path.join(args.slivar_data_dir, "imprinting.txt")), "Gene")
    pseudoautosomal = index_first_by(load_reference_table(os.path.join(args.slivar_data_dir, "pseudoautosomal.txt")), "Ensembl_gene_id")
    hgmd_by_variant, hgmd_genes = load_hgmd(args.hgmd)

    rows = []
    with pysam.VariantFile(args.vcf) as vcf:
        samples = list(vcf.header.samples)
        # Map original VCF sample IDs to report column suffixes; pysam still uses the original IDs.
        # But reporting without hyphens.
        sample_headers = {sample: report_sample_name(sample) for sample in samples}
        # DN and DQ are optional FORMAT fields; include those columns only when present.
        include_denovo = "DN" in vcf.header.formats
        include_denovo_quality = "DQ" in vcf.header.formats
        columns = make_columns(args.mode, samples, include_denovo, include_denovo_quality, args.profile)
        csq_fields = get_csq_fields(vcf)
        if not csq_fields:
            raise SystemExit("VCF header does not contain a usable INFO/CSQ field")

        # postfilter.py has already decomposed and deduplicated these variant records.
        for record in vcf:
            if not record.alts:
                continue
            csq_records = parse_csq_annotations(record, csq_fields, start_index=1)
            # Choose the final report CSQ by prioritizing non-pseudogene, impactful variations,
            # protein-coding, variation order, named-gene, MANE/canonical annotations.
            primary = select_primary_csq(csq_records, order_map)

            ref = record.ref
            alt = record.alts[0]
            position = f"{record.chrom}:{record.pos}"

            row = {column: "" for column in columns}
            row["Position"] = position
            row["UCSC_Link"] = hyperlink_ucsc(position)
            row["GNOMAD_Link"] = hyperlink_gnomad(record.chrom, record.pos, ref, alt)
            row["Ref"] = ref
            row["Alt"] = alt

            # These lists become the combined gts and Trio_coverage fields after the sample columns.
            sample_gt_strings = []
            sample_depths = []
            for sample in samples:
                # Use the original sample ID for pysam, then the report-safe suffix for CSV columns.
                sample_data = record.samples[sample]
                sample_header = sample_headers[sample]
                sample_gt_strings.append(format_genotype(sample_data, ref, record.alts or []))
                sample_depths.append(get_sample_depth(sample_data))
                alt_depth = get_alt_depth(sample_data)
                sample_gq = get_sample_gq(sample_data)
                row[f"Zygosity.{sample_header}"] = classify_zygosity(sample_data, ref, record.alts or [], record.chrom)
                row[f"Alt_depths.{sample_header}"] = "" if alt_depth is None else str(alt_depth)
                # Missing genotype quality is reported as -1 to match the CRE reports.
                if has_value(sample_gq):
                    row[f"gt_quals.{sample_header}"] = sample_gq
                else:
                    row[f"gt_quals.{sample_header}"] = "-1"
                # DN/DQ are optional FORMAT fields, so their columns are only filled when present.
                if include_denovo:
                    row[f"denovo.{sample_header}"] = value_as_text(sample_data.get("DN"))
                if include_denovo_quality:
                    row[f"denovo_quality.{sample_header}"] = value_as_text(sample_data.get("DQ"))

            row["gts"] = ",".join(sample_gt_strings)
            if args.profile == "pacbio":
                row["PS"] = format_phase_sets(get_info_value(record, "PS"))
            row["Depth"] = value_as_text(get_info_value(record, "DP"))
            row["Quality"] = "" if record.qual is None else str(record.qual)
            row["Trio_coverage"] = "_".join(depth if has_value(depth) else "0" for depth in sample_depths)
            row["Clinvar"] = clinvar_text(record) or "."
            hgmd_match = hgmd_by_variant.get(f"{position}-{ref}-{alt}", {})
            row["HGMD_id"] = hgmd_match.get("HGMD_id", "NA")
            # HGMD_gene is filled later once the primary and overlapping genes are known.
            row["HGMD_gene"] = "NA"
            row["HGMD_tag"] = hgmd_match.get("HGMD_tag", "NA")
            row["HGMD_ref"] = hgmd_match.get("HGMD_ref", "NA")

            # Copy direct INFO annotations into report columns, with special handling
            # for max-per-variant scores, binary flags, and aligned GreenDB lists.
            for field, source in [
                ("Gnomad_af_grpmax", "gnomad_af_grpmax"),
                ("Gnomad_af", "gnomad_af"), ("Gnomad_ac", "gnomad_ac"),
                ("Gnomad_hom", "gnomad_hom"), ("Gnomad_male_ac", "gnomad_male_ac"),
                ("Gnomad_fafmax_faf95_max", "gnomad_fafmax_faf95_max"),
                ("Gnomad_filter", "gnomad_filter"),
                ("Regeneron_exome_AF", "regeneron_exome_AF"), ("Regeneron_exome_AC", "regeneron_exome_AC"),
                ("thousandG_AF", "thousandG_AF"), ("thousandG_AC", "thousandG_AC"),
                ("thousandG_nhomalt", "thousandG_nhomalt"),
                ("GSO_AF", "GSO_AF"), ("GSO_AC", "GSO_AC"),
                ("GSO_nhomalt", "GSO_nhomalt"), ("GSO_hemi", "GSO_hemi"),
                ("CoLoRSdb_AF", "CoLoRSdb_AF"), ("CoLoRSdb_AC", "CoLoRSdb_AC"),
                ("CoLoRSdb_AC_Hemi", "CoLoRSdb_AC_Hemi"), ("CoLoRSdb_nhomalt", "CoLoRSdb_nhomalt"),
                ("TG_LRWGS_AC", "tg_lrwgs_ac"), ("TG_LRWGS_samples", "tg_lrwgs_samples"),
                ("TG_LRWGS_hom", "tg_lrwgs_hom"), ("Dark_genes", "Dark_genes"),
                ("rsIDs", "rs_ids"),
                ("Cadd_score", "CADD_phred"), ("Vest4_score", "Vest4_score"),
                ("Revel_score", "REVEL_score"), ("Gerp_score", "Gerp_score"),
                ("AlphaMissense", "AlphaMissense"),
                ("ncER_score", "ncER"), ("ReMM_score", "ReMM"),
                ("LINSIGHT_score", "LinSight_Score"),
                ("TF_binding_sites", "tf_binding_sites"),
                ("GreenDB_variant_type", "GreenDB_variant_type"),
                ("GreenDB_closest_gene", "GreenDB_closest_gene"), ("GreenDB_controlled_gene", "GreenDB_controlled_gene"),
                ("CTCF_binding_site", "CTCF_binding_site"), ("DNaseI_hypersensitive_site", "DNaseI_hypersensitive_site"),
                ("ENH_cellline_tissue", "ENH_cellline_tissue"),
                ("UCE_100bp", "UCE_100bp"), ("UCE_200bp", "UCE_200bp"),
                ("Old_multiallelic", "OLD_MULTIALLELIC"),
            ]:
                value = get_info_value(record, source)
                if field == "Vest4_score":
                    # VEST4 can contain several comma-separated scores; report the highest one.
                    set_if_column_present(row, field, max_numeric_csv(value))
                elif field in FLAG_BINARY_FIELDS:
                    # Presence/absence INFO flags become explicit 1/0 report values.
                    set_if_column_present(row, field, "1" if has_value(value) else "0")
                elif field in ALIGNED_LIST_FIELDS:
                    # GreenDB region columns must preserve empty slots across aligned lists.
                    set_if_column_present(row, field, as_aligned_text(value))
                else:
                    set_if_column_present(row, field, value_as_text(value))

            # phyloP100way is a CSQ/transcript value; prefer the selected CSQ, then any CSQ.
            phylop100way = ""
            if primary is not None:
                phylop100way = primary.get("phyloP100way", "")
            if not has_value(phylop100way):
                for csq in csq_records:
                    phylop100way = csq.get("phyloP100way", "")
                    if has_value(phylop100way):
                        break
            row["phylop100way"] = value_as_text(phylop100way)

            # Use the missing token for absent filters/scores before final normalization.
            if "Gnomad_filter" in row:
                row["Gnomad_filter"] = row["Gnomad_filter"] or "None"
            if "Cadd_score" in row:
                row["Cadd_score"] = row["Cadd_score"] or "None"
            for field in ["Vest4_score", "Revel_score", "Gerp_score", "AlphaMissense", "ncER_score", "ReMM_score", "LINSIGHT_score"]:
                if field in row:
                    row[field] = row.get(field, "") or "None"

            spliceai_impact, spliceai_score = parse_spliceai(get_info_value(record, "spliceai_score"))
            set_if_column_present(row, "SpliceAI_impact", spliceai_impact)
            set_if_column_present(row, "SpliceAI_score", spliceai_score)
            set_if_column_present(row, "promoterAI_score", parse_promoterai(get_info_value(record, "promoterAI")))
            if "Noncoding_path_pred" in row:
                row["Noncoding_path_pred"] = noncoding_pred_fraction(
                    row.get("Cadd_score", ""), row.get("ncER_score", ""), row.get("ReMM_score", ""), row.get("LINSIGHT_score", "")
                )

            # Build *_all columns from every CSQ annotation before choosing what each report field needs.
            all_gene_symbols = []
            all_gene_ids = []
            # Some CSQs have an Ensembl gene ID but no gene symbol. The loop below
            # keeps those IDs without adding a blank/made-up Gene_all value.
            unnamed_gene_ids = []
            all_ensembl_transcripts = []
            variation_all_values = []
            sift_all_values = []
            polyphen_all_values = []
            biotype_all_values = []
            impact_all_values = []

            # Split CSQ genes into named entries for Gene_all and ID-only entries for Ensembl_gene_id_all.
            for csq in csq_records:
                symbol = get_gene_symbol(csq)
                if has_value(symbol):
                    # Named CSQs contribute a display symbol and, when available, a matching gene ID.
                    all_gene_symbols.append(symbol)
                    gene_id = select_ensembl_gene_id(csq_records, csq)
                    if has_value(gene_id):
                        all_gene_ids.append(gene_id)
                else:
                    # Symbol-less CSQs only contribute an Ensembl gene ID.
                    ensembl_gene_id = get_csq_ensembl_gene_id(csq)
                    if has_value(ensembl_gene_id):
                        unnamed_gene_ids.append(ensembl_gene_id)

                feature = csq.get("Feature", "")
                if feature.startswith("ENST"):
                    all_ensembl_transcripts.append(feature)

                # Keep raw per-CSQ values here; join_unique_values below applies de-dup/order.
                variation_all_values.append(csq.get("Consequence", ""))
                sift_all_values.append(parse_vep_score(csq.get("SIFT", "")))
                polyphen_all_values.append(parse_vep_score(csq.get("PolyPhen", "")))
                biotype_all_values.append(csq.get("BIOTYPE", ""))
                impact_all_values.append(csq.get("IMPACT", ""))

            # Collapse per-CSQ lists into the comma-delimited *_all report strings.
            row["Gene_all"] = join_unique_values(all_gene_symbols)
            row["Ensembl_gene_id_all"] = join_unique_values(all_gene_ids + unnamed_gene_ids)
            row["Ensembl_transcript_id_all"] = join_unique_values(all_ensembl_transcripts)
            row["Variation_all"] = join_unique_values(variation_all_values)
            row["Info_all"] = join_unique_values(build_info_item(csq) for csq in csq_records)
            row["Refseq_change_all"] = refseq_change_all(csq_records)
            set_if_column_present(row, "Sift_score_all", join_unique_values(sift_all_values) or "None")
            set_if_column_present(row, "Polyphen_score_all", join_unique_values(polyphen_all_values) or "None")
            row["CSQ_biotype_all"] = join_unique_values(biotype_all_values)
            row["CSQ_impact_all"] = join_unique_values(impact_all_values)

            gene_desc_all = []
            omim_pheno_all = []
            omim_inh_all = []
            orphanet_all = []
            # Gene-level reference joins use Gene_all symbols, with ENSG used for
            # gene descriptions and Orphanet when available.
            for symbol in row["Gene_all"].split(","):
                if not has_value(symbol):
                    continue
                # Gene_all is symbol text; these reference tables are keyed by Ensembl gene ID.
                ensg = find_ensembl_gene_id(csq_records, symbol)
                if has_value(ensg):
                    gene_desc_all.append(gene_descriptions.get(ensg, {}).get("Gene_description", ""))
                    orphanet_all.append(lookup_orphanet(ensg, orphanet))
                phen, inh = join_omim(symbol, omim)
                omim_pheno_all.append(phen)
                omim_inh_all.append(inh)
            row["Gene_description_all"] = join_unique_values(gene_desc_all)
            row["omim_phenotype_all"] = join_unique_values(omim_pheno_all)
            row["omim_inheritance_all"] = join_unique_values(omim_inh_all)
            row["Orphanet_all"] = join_unique_values(orphanet_all)

            # Constraint_all summarizes every Ensembl transcript found in the CSQs.
            row["Constraint_all"] = constraint_all_summary(row["Ensembl_transcript_id_all"].split(","), constraint)

            # Fill single-transcript report fields from the selected primary CSQ.
            if primary is not None:
                primary_gene = get_gene_symbol(primary)
                primary_ensg = select_ensembl_gene_id(csq_records, primary)
                primary_feature = primary.get("Feature", "")
                # The constraint table is keyed by Ensembl transcript IDs.
                primary_tx = primary_feature if primary_feature.startswith("ENST") else ""
                row["Gene"] = primary_gene
                set_if_column_present(row, "Ensembl_gene_id", primary_ensg)
                row["Variation"] = select_consequence(primary, order_map)
                row["Info"] = build_info_item(primary) or "NA"
                # Info follows the selected primary CSQ; Refseq_change separately chooses
                # the best RefSeq annotation on that gene, similar in intent to CRE's NM_ scan.
                row["Refseq_change"] = best_refseq_change_for_gene(csq_records, primary_gene, order_map)
                row["Ensembl_transcript_id"] = primary_tx
                set_if_column_present(row, "AA_position", replace_slashes(primary.get("Protein_position", "")))
                set_if_column_present(row, "Exon", replace_slashes(primary.get("EXON", "")))
                set_if_column_present(row, "Protein_domains", primary.get("DOMAINS", ""))
                set_if_column_present(row, "Sift_score", parse_vep_score(primary.get("SIFT", "")) or "None")
                set_if_column_present(row, "Polyphen_score", parse_vep_score(primary.get("PolyPhen", "")) or "None")
                row["CSQ_biotype"] = primary.get("BIOTYPE", "")
                row["CSQ_impact"] = primary.get("IMPACT", "")
                # Gene-symbol lookups try the primary gene first, then other Gene_all entries.
                overlapping_genes = genes_primary_first(primary_gene, row["Gene_all"])
                # Imprinting and HGMD can belong to an overlapping gene rather than the selected CSQ gene.
                row["Imprinting_status"], row["Imprinting_expressed_allele"] = lookup_imprinting(overlapping_genes, imprinting)
                # Pseudoautosomal is keyed by ENSG, so use the CSQ gene IDs gathered above.
                row["Pseudoautosomal"] = lookup_pseudoautosomal(all_gene_ids, pseudoautosomal)
                row["HGMD_gene"] = first_hgmd_gene(overlapping_genes, hgmd_genes)
                # Constraint columns are transcript-level and only use the selected Ensembl transcript.
                for field, value in primary_constraint_values(primary_tx, constraint).items():
                    row[field] = value
            else:
                row["Info"] = "NA"
                row["Refseq_change"] = "NA"
                row["HGMD_gene"] = "NA"

            # Normalize missing-value tokens after every row field has had a chance to be set.
            for field in row:
                row[field] = normalize_report_value(field, row[field])

            rows.append(row)

    # Burden is the per-sample count of unique Het/Hom variants for each gene.
    # Key by variant+gene so duplicate rows do not inflate the count.
    burden_seen = {sample: set() for sample in samples}
    burden_by_sample = {sample: Counter() for sample in samples}
    for row in rows:
        gene = row["Gene"]
        if not has_value(gene):
            continue
        variant_gene_key = (row["Position"], row["Ref"], row["Alt"], gene)
        for sample in samples:
            # The key prevents duplicate rows for the same variant/gene from inflating burden.
            if row[f"Zygosity.{sample_headers[sample]}"] in {"Het", "Hom"} and variant_gene_key not in burden_seen[sample]:
                burden_seen[sample].add(variant_gene_key)
                burden_by_sample[sample][gene] += 1
    for row in rows:
        for sample in samples:
            row[f"Burden.{sample_headers[sample]}"] = str(burden_by_sample[sample].get(row["Gene"], 0))

    columns = drop_empty_optional_columns(columns, rows, ("denovo.", "denovo_quality."))

    # Match the previous grouped-record order for variants sharing one position.
    rows.sort(key=lambda row: (row["Position"], row["Ref"], row["Alt"]))
    with open(args.out_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
