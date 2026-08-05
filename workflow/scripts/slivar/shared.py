#!/usr/bin/env python3
# add docstring for functions

import os

MISSING = {"", ".", "None", "NA"}

# Load the curated RNA disease-gene list (one symbol per line; "#" comments).
def _load_disease_rna_genes(path):
    genes = set()
    with open(path) as handle:
        for line in handle:
            symbol = line.strip()
            if not symbol or symbol.startswith("#"):
                continue
            genes.add(symbol)
    return genes


# Small non-coding disease genes VEP annotates under an overlapping protein-coding
# gene. Membership promotes a CSQ to the primary report gene.
DISEASE_RNA_GENES = _load_disease_rna_genes(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "disease-rna-genes.txt")
)


# Check whether a scalar, boolean, or tuple INFO value contains real content.
def has_value(value):
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, tuple):
        return any(has_value(item) for item in value)
    return str(value) not in MISSING


def get_info_value(record, field):
    try:
        return record.info.get(field)
    except (KeyError, ValueError):
        return None


def value_as_text(value):
    if value is None:
        return ""
    if isinstance(value, tuple):
        return ",".join(str(item) for item in value if has_value(item))
    text = str(value)
    if text in MISSING:
        return ""
    return text


# Convert scalar INFO values, or the first value of a tuple INFO field, to float.
def value_as_float(value):
    if value is None:
        return None
    if isinstance(value, tuple):
        if len(value) == 0:
            return None
        value = value[0]
    if str(value) in MISSING:
        return None
    try:
        return float(value)
    except Exception:
        return None

# Remove "_variant" from VCF consequence names to match the order-file terms.
def normalize_consequence_term(term):
    text = str(term).strip()
    if not text:
        return ""
    if text.endswith("_variant"):
        return text[: -len("_variant")]
    return text

# Load consequence terms as a normalized rank map.
def load_consequence_order(path):
    order = {}
    rank = 0
    with open(path) as handle:
        for line in handle:
            term = normalize_consequence_term(line.strip())
            # Ignore blank lines and lines beginning with "#".
            if not term or term.startswith("#"):
                continue
            order[term] = rank
            rank += 1
    if not order:
        raise SystemExit(f"No consequence terms found in {path}")
    return order


# Extract pipe-delimited VEP CSQ field names from the INFO/CSQ header.
def get_csq_fields(vcf):
    if "CSQ" not in vcf.header.info:
        return []
    description = vcf.header.info["CSQ"].description
    marker = "Format:"
    if marker not in description:
        return []
    return description.split(marker, 1)[1].strip().strip('"').split("|")


# Parse a record's INFO/CSQ values into one dict per VEP annotation.
# csq_fields comes from the CSQ header; each raw annotation is pipe-delimited
# in the same order. Missing trailing fields are padded, and _csq_index
# preserves input order for downstream tie-breaking.
def parse_csq_annotations(record, csq_fields, start_index=0):
    if not csq_fields:
        return []
    raw = get_info_value(record, "CSQ")
    if raw is None:
        return []
    # pysam usually exposes multiple comma-separated CSQ annotations as a tuple.
    if isinstance(raw, tuple):
        entries = raw
    else:
        entries = [raw]
    annotations = []
    for index, entry in enumerate(entries, start=start_index):
        # Match each pipe-delimited value to the corresponding CSQ header field.
        values = str(entry).split("|")
        if len(values) < len(csq_fields):
            # Pad truncated annotations so later fields still map to empty strings.
            values.extend([""] * (len(csq_fields) - len(values)))
        # VEP CSQ values follow the header order, so zip pairs names with values.
        annotation = dict(zip(csq_fields, values))
        annotation["_csq_index"] = index
        annotations.append(annotation)
    return annotations


# VEP can join multiple consequences with "&"; return the trimmed terms.
def get_consequence_terms(csq):
    consequence = csq.get("Consequence", "")
    if not has_value(consequence):
        return []
    terms = []
    for term in consequence.split("&"):
        term = term.strip()
        if term:
            terms.append(term)
    return terms


# Normalize each consequence term before comparing against the order file.
def _normalized_consequence_terms(csq):
    normalized_terms = []
    for term in get_consequence_terms(csq):
        normalized_term = normalize_consequence_term(term)
        if normalized_term:
            normalized_terms.append(normalized_term)
    return normalized_terms


# Return the best-ranked consequence index for a CSQ; unknown terms sort last.
def rank_consequence(csq, consequence_order):
    # One CSQ can have multiple "&"-joined consequences. Convert each known
    # consequence to its order-file rank, then keep the lowest rank.
    ranks = []
    # Get "Consequence" from the CSQ annotation, normalize the terms, and rank them.
    for term in _normalized_consequence_terms(csq):
        if term in consequence_order:
            ranks.append(consequence_order[term])
    if ranks:
        return min(ranks)
    # Unknown-only consequences should sort after every ranked consequence.
    return 10**9


def select_consequence(csq, consequence_order):
    terms = get_consequence_terms(csq)
    if not terms:
        return ""
    return sorted(
        terms,
        key=lambda term: (
            consequence_order.get(normalize_consequence_term(term), 10**9),
            term,
        ),
    )[0]


# Return a display gene symbol for a CSQ. SYMBOL/HGNC are preferred; Gene is
# used only when it is already a symbol rather than an Ensembl gene ID.
def get_gene_symbol(csq):
    if csq is None:
        return ""
    for field in ("SYMBOL", "HGNC"):
        value = csq.get(field, "")
        if has_value(value):
            return value
    gene = csq.get("Gene", "")
    if has_value(gene) and not str(gene).startswith("ENSG"):
        return gene
    return ""


def get_csq_ensembl_gene_id(csq):
    gene = csq.get("Gene", "")
    if gene.startswith("ENSG"):
        return gene
    return ""


# Translate a report gene symbol back to the first matching Ensembl gene ID in the CSQs.
def find_ensembl_gene_id(csq_annotations, symbol):
    candidates = []
    for csq in csq_annotations:
        if get_gene_symbol(csq) != symbol:
            continue
        gene = get_csq_ensembl_gene_id(csq)
        if gene:
            # Keep candidates in CSQ/input order so duplicate symbols resolve consistently.
            candidates.append(gene)
    if candidates:
        return candidates[0]
    return ""


def select_ensembl_gene_id(csq_annotations, primary_csq):
    if primary_csq is None:
        return ""
    symbol = get_gene_symbol(primary_csq)
    ensembl_gene_id = find_ensembl_gene_id(csq_annotations, symbol)
    if has_value(ensembl_gene_id):
        return ensembl_gene_id
    gene = primary_csq.get("Gene", "")
    if has_value(gene):
        return gene
    return ""

# Return which variations are above a defined cutoff in the order.
def _above_cutoff_rank(csq, consequence_order, cutoff):
    cutoff_rank = consequence_order.get(cutoff)
    if cutoff_rank is None:
        return 1
    if rank_consequence(csq, consequence_order) < cutoff_rank:
        return 0
    return 1


# Treat placeholder-style symbols as weaker than named genes for report display.
def _is_weak_gene_symbol(symbol):
    if not has_value(symbol):
        return False
    text = str(symbol).upper()
    if text.startswith(("ENSG", "LOC", "LINC", "MIR", "RN7", "RNU", "SNORD", "SNORA", "SCARNA", "Y_RNA")):
        return True
    for prefix in ("AC", "AL", "AP"):
        if text.startswith(prefix) and len(text) > len(prefix) and text[len(prefix)].isdigit():
            return True
    return text.startswith("RP") and len(text) > 2 and (text[2].isdigit() or text[2] == "-")


# Return the tuple used to choose the final report CSQ for display fields.
def primary_csq_sort_key(csq, consequence_order):
    biotype = str(csq.get("BIOTYPE", "")).lower()
    symbol = get_gene_symbol(csq)

    # Promote curated RNA disease genes ahead of every other tie-breaker. VEP would
    # otherwise pick the overlapping protein-coding gene as primary and bury the RNA
    # gene. Display-only: this does not affect which variants are kept.
    if has_value(symbol) and symbol in DISEASE_RNA_GENES:
        curated_rna_rank = 0
    else:
        curated_rna_rank = 1

    # Prefer non-pseudogene annotations before other transcript tie-breakers.
    if "pseudogene" in biotype:
        pseudogene_rank = 1
    else:
        pseudogene_rank = 0
    # IMPACTFUL_CUTOFF is a sentinel in the consequence order file.
    impactful_rank = _above_cutoff_rank(csq, consequence_order, "IMPACTFUL_CUTOFF")

    # Keep protein-coding annotations ahead of related decay classes, then other biotypes.
    if biotype == "protein_coding":
        protein_coding_rank = 0
    elif biotype.startswith("protein_coding_") or biotype in {"nonsense_mediated_decay", "non_stop_decay"}:
        protein_coding_rank = 1
    else:
        protein_coding_rank = 2

    # Within the same broad class, use the order-file consequence severity.
    consequence_rank = rank_consequence(csq, consequence_order)
    # Prefer named gene symbols, then placeholder-style symbols, then missing symbols.
    if not has_value(symbol):
        gene_rank = 2
    elif _is_weak_gene_symbol(symbol):
        gene_rank = 1
    else:
        gene_rank = 0

    # MANE and canonical flags are transcript-quality tie-breakers.
    if has_value(csq.get("MANE_SELECT", "")):
        mane_rank = 0
    elif has_value(csq.get("MANE_PLUS_CLINICAL", "")):
        mane_rank = 1
    else:
        mane_rank = 2

    if csq.get("CANONICAL", "") == "YES":
        canonical_rank = 0
    else:
        canonical_rank = 1

    if biotype == "processed_transcript":
        transcript_rank = 0
    elif biotype in {"nonsense_mediated_decay", "non_stop_decay", "retained_intron"} or biotype.endswith("_transcript"):
        transcript_rank = 1
    else:
        transcript_rank = 2

    # GENIC_CUTOFF is another order-file sentinel; it breaks late ties toward genic CSQs.
    genic_rank = _above_cutoff_rank(csq, consequence_order, "GENIC_CUTOFF")

    # min() compares tuple entries from left to right, so lower ranks win in
    # this order: curated RNA gene, non-pseudogene, impactful, protein-coding,
    # consequence, gene symbol quality, MANE, canonical, transcript class, genic,
    # transcript ID, input order.
    return (
        curated_rna_rank,
        pseudogene_rank,
        impactful_rank,
        protein_coding_rank,
        consequence_rank,
        gene_rank,
        mane_rank,
        canonical_rank,
        transcript_rank,
        genic_rank,
        csq.get("Feature", ""),
        csq.get("_csq_index", 0),
    )


# Select the CSQ used for final report fields.
def select_primary_csq(csq_annotations, consequence_order):
    if not csq_annotations:
        return None
    return min(
        csq_annotations,
        key=lambda csq: primary_csq_sort_key(csq, consequence_order),
    )


# Convert pysam FORMAT/GT allele indexes into the allele strings used in reports.
def format_genotype(sample_data, ref, alts):
    # pysam stores GT as allele indexes: 0 is REF, 1+ are entries in ALTS, and None is missing.
    gt = sample_data.get("GT")
    if gt is None:
        return "./."
    # If every allele is missing, report the same missing genotype string CRE used.
    if all(allele_index is None for allele_index in gt):
        return "./."

    alleles = []
    for allele_index in gt:
        if allele_index is None:
            alleles.append(".")
        elif allele_index == 0:
            alleles.append(ref)
        else:
            try:
                # ALT indexes are one-based in GT, so subtract one for the Python list.
                alleles.append(alts[allele_index - 1])
            except Exception:
                alleles.append(".")
    # Keep the VCF separator for diploid calls, including phased half-calls.
    if sample_data.phased:
        separator = "|"
    else:
        separator = "/"
    return separator.join(alleles)


# Return the strongest ALT support from FORMAT/AD.
def get_alt_depth(sample_data):
    ad = sample_data.get("AD")
    if ad is None:
        return None

    # AD is REF depth followed by one depth per ALT allele; skip the REF depth.
    if isinstance(ad, str):
        alt_depths = ad.split(",")[1:]
    else:
        try:
            alt_depths = list(ad)[1:]
        except Exception:
            return None

    parsed = []
    for depth in alt_depths:
        # Missing or malformed ALT depths are ignored instead of counted as zero.
        if depth is None or str(depth) in MISSING:
            continue
        try:
            parsed.append(int(depth))
        except Exception:
            continue
    if parsed:
        # Multiallelic records can have several ALT depths; keep the largest one.
        return max(parsed)
    return None


def get_sample_depth(sample_data):
    depth = sample_data.get("DP")
    if depth is None:
        return ""
    return str(depth)


def get_sample_gq(sample_data):
    gq = sample_data.get("GQ")
    if gq is None:
        return ""
    return str(gq)
