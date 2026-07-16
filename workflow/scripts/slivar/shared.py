#!/usr/bin/env python3


MISSING = {"", ".", "None", "NA"}


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
    return "" if text in MISSING else text


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


# Consequence ordering and CSQ selection shared by report and CH generation.
def normalize_consequence_term(term):
    text = str(term).strip()
    if not text:
        return ""
    if text.endswith("_variant"):
        return text[: -len("_variant")]
    return text


def load_consequence_order(path):
    order = {}
    rank = 0
    with open(path) as handle:
        for line in handle:
            term = normalize_consequence_term(line.strip())
            if not term or term.startswith("#"):
                continue
            order[term] = rank
            rank += 1
    if not order:
        raise SystemExit(f"No consequence terms found in {path}")
    return order


def get_csq_fields(vcf):
    if "CSQ" not in vcf.header.info:
        return []
    description = vcf.header.info["CSQ"].description
    marker = "Format:"
    if marker not in description:
        return []
    return description.split(marker, 1)[1].strip().strip('"').split("|")


def parse_csq_annotations(record, csq_fields, start_index=0):
    if not csq_fields:
        return []
    raw = get_info_value(record, "CSQ")
    if raw is None:
        return []
    entries = raw if isinstance(raw, tuple) else [raw]
    annotations = []
    for index, entry in enumerate(entries, start=start_index):
        values = str(entry).split("|")
        if len(values) < len(csq_fields):
            values.extend([""] * (len(csq_fields) - len(values)))
        annotation = dict(zip(csq_fields, values))
        annotation["_csq_index"] = index
        annotations.append(annotation)
    return annotations


def get_consequence_terms(csq):
    consequence = csq.get("Consequence", "")
    if not has_value(consequence):
        return []
    return [term.strip() for term in consequence.split("&") if term.strip()]


def _normalized_consequence_terms(csq):
    return [
        normalize_consequence_term(term)
        for term in get_consequence_terms(csq)
        if normalize_consequence_term(term)
    ]


def rank_consequence(csq, consequence_order):
    ranks = [
        consequence_order[term]
        for term in _normalized_consequence_terms(csq)
        if term in consequence_order
    ]
    return min(ranks) if ranks else 10**9


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
    return gene if gene.startswith("ENSG") else ""


def find_ensembl_gene_id(csq_annotations, symbol):
    candidates = []
    for csq in csq_annotations:
        if get_gene_symbol(csq) != symbol:
            continue
        gene = get_csq_ensembl_gene_id(csq)
        if gene:
            candidates.append(gene)
    return candidates[0] if candidates else ""


def select_ensembl_gene_id(csq_annotations, primary_csq):
    if primary_csq is None:
        return ""
    symbol = get_gene_symbol(primary_csq)
    ensembl_gene_id = find_ensembl_gene_id(csq_annotations, symbol)
    if has_value(ensembl_gene_id):
        return ensembl_gene_id
    gene = primary_csq.get("Gene", "")
    return gene if has_value(gene) else ""


def _above_cutoff_rank(csq, consequence_order, cutoff):
    cutoff_rank = consequence_order.get(cutoff)
    if cutoff_rank is None:
        return 1
    return 0 if rank_consequence(csq, consequence_order) < cutoff_rank else 1


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


def primary_csq_sort_key(csq, consequence_order):
    biotype = str(csq.get("BIOTYPE", "")).lower()
    symbol = get_gene_symbol(csq)

    pseudogene_rank = 1 if "pseudogene" in biotype else 0
    impactful_rank = _above_cutoff_rank(csq, consequence_order, "IMPACTFUL_CUTOFF")

    if biotype == "protein_coding":
        protein_coding_rank = 0
    elif biotype.startswith("protein_coding_") or biotype in {"nonsense_mediated_decay", "non_stop_decay"}:
        protein_coding_rank = 1
    else:
        protein_coding_rank = 2

    consequence_rank = rank_consequence(csq, consequence_order)
    gene_rank = 2 if not has_value(symbol) else (1 if _is_weak_gene_symbol(symbol) else 0)

    if has_value(csq.get("MANE_SELECT", "")):
        mane_rank = 0
    elif has_value(csq.get("MANE_PLUS_CLINICAL", "")):
        mane_rank = 1
    else:
        mane_rank = 2

    canonical_rank = 0 if csq.get("CANONICAL", "") == "YES" else 1

    if biotype == "processed_transcript":
        transcript_rank = 0
    elif biotype in {"nonsense_mediated_decay", "non_stop_decay", "retained_intron"} or biotype.endswith("_transcript"):
        transcript_rank = 1
    else:
        transcript_rank = 2

    genic_rank = _above_cutoff_rank(csq, consequence_order, "GENIC_CUTOFF")

    return (
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


def select_primary_csq(csq_annotations, consequence_order):
    if not csq_annotations:
        return None
    return min(
        csq_annotations,
        key=lambda csq: primary_csq_sort_key(csq, consequence_order),
    )


# Genotype values that must agree between report and CH output.
def format_genotype(sample_data, ref, alts):
    gt = sample_data.get("GT")
    if gt is None:
        return "./."
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
                alleles.append(alts[allele_index - 1])
            except Exception:
                alleles.append(".")
    # Keep the VCF separator for diploid calls, including phased half-calls.
    return ("|" if sample_data.phased else "/").join(alleles)


def get_alt_depth(sample_data):
    ad = sample_data.get("AD")
    if ad is None:
        return None

    if isinstance(ad, str):
        alt_depths = ad.split(",")[1:]
    else:
        try:
            alt_depths = list(ad)[1:]
        except Exception:
            return None

    parsed = []
    for depth in alt_depths:
        if depth is None or str(depth) in MISSING:
            continue
        try:
            parsed.append(int(depth))
        except Exception:
            continue
    return max(parsed) if parsed else None


def get_sample_depth(sample_data):
    depth = sample_data.get("DP")
    return "" if depth is None else str(depth)


def get_sample_gq(sample_data):
    gq = sample_data.get("GQ")
    return "" if gq is None else str(gq)
