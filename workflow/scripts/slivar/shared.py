#!/usr/bin/env python3

import csv
from collections import defaultdict


MISSING = {"", ".", "None", "NA"}


def present(value):
    if value is None:
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, tuple):
        return any(present(v) for v in value)
    return str(value) not in MISSING


def info(record, field):
    try:
        return record.info.get(field)
    except (KeyError, ValueError):
        return None


def as_text(value):
    if value is None:
        return ""
    if isinstance(value, tuple):
        return ",".join(str(v) for v in value if present(v))
    text = str(value)
    return "" if text in MISSING else text


def as_float(value):
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


def uniq_join(values, sep=","):
    seen = set()
    out = []
    for value in values:
        if value is None:
            continue
        items = value if isinstance(value, (tuple, list)) else [value]
        for item in items:
            text = str(item)
            if text in MISSING or text in seen:
                continue
            seen.add(text)
            out.append(text)
    return sep.join(out)


def clean_slash(value):
    if value is None:
        return ""
    return str(value).replace("/", "_")


def load_table(path, delimiter=None):
    if not path:
        return []
    if delimiter is None:
        delimiter = "," if path.endswith(".csv") else "\t"
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def normalize_consequence_term(term):
    text = str(term).strip()
    if not text:
        return ""
    if text.endswith("_variant"):
        return text[: -len("_variant")]
    return text


def load_impact_order(path):
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


def index_first(rows, key_column):
    index = {}
    for row in rows:
        key = row.get(key_column, "")
        if key and key not in index:
            index[key] = row
    return index


def index_many(rows, key_column):
    index = defaultdict(list)
    for row in rows:
        key = row.get(key_column, "")
        if key:
            index[key].append(row)
    return index


def index_constraint(rows):
    index = {}
    for row in rows:
        transcript = row.get("Ensembl_transcript_id", "")
        if not transcript:
            continue
        index[transcript] = row
        transcript_no_version = transcript.split(".", 1)[0]
        if transcript_no_version not in index:
            index[transcript_no_version] = row
    return index


def parse_csq_header(vcf):
    if "CSQ" not in vcf.header.info:
        return []
    description = vcf.header.info["CSQ"].description
    marker = "Format:"
    if marker not in description:
        return []
    return description.split(marker, 1)[1].strip().strip('"').split("|")


def parse_csq_records(record, csq_fields, start_index=0):
    if not csq_fields:
        return []
    raw = info(record, "CSQ")
    if raw is None:
        return []
    entries = raw if isinstance(raw, tuple) else [raw]
    parsed = []
    for idx, entry in enumerate(entries, start=start_index):
        values = str(entry).split("|")
        if len(values) < len(csq_fields):
            values.extend([""] * (len(csq_fields) - len(values)))
        csq = dict(zip(csq_fields, values))
        csq["_csq_index"] = idx
        parsed.append(csq)
    return parsed


def consequence_terms(csq):
    consequence = csq.get("Consequence", "")
    if not present(consequence):
        return []
    return [term.strip() for term in consequence.split("&") if term.strip()]


def normalized_consequence_terms(csq):
    return [normalize_consequence_term(term) for term in consequence_terms(csq) if normalize_consequence_term(term)]


def gene_symbol(csq):
    if csq is None:
        return ""
    for field in ("SYMBOL", "HGNC"):
        value = csq.get(field, "")
        if present(value):
            return value
    gene = csq.get("Gene", "")
    if present(gene) and not str(gene).startswith("ENSG"):
        return gene
    return ""


def cre_report_gene(csq):
    # High-impact filtering should use the same value shown in the report Gene column.
    return gene_symbol(csq)


def consequence_rank(csq, order_map):
    ranks = [order_map[term] for term in normalized_consequence_terms(csq) if term in order_map]
    return min(ranks) if ranks else 10**9


def consequence_display(csq, order_map):
    terms = consequence_terms(csq)
    if not terms:
        return ""
    return sorted(terms, key=lambda term: (order_map.get(normalize_consequence_term(term), 10**9), term))[0]


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


def gene_id_with_fallback_for_csq(csq_records, primary):
    if primary is None:
        return ""
    symbol = gene_symbol(primary)
    ensg = best_ensembl_gene_id_for_symbol(csq_records, symbol)
    if present(ensg):
        return ensg
    gene = primary.get("Gene", "")
    return gene if present(gene) else ""


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


def pseudogene_rank(csq):
    biotype = str(csq.get("BIOTYPE", "")).lower()
    return 1 if "pseudogene" in biotype else 0


def protein_coding_biotype_rank(csq):
    biotype = str(csq.get("BIOTYPE", "")).lower()
    if biotype == "protein_coding":
        return 0
    if biotype.startswith("protein_coding_") or biotype in {"nonsense_mediated_decay", "non_stop_decay"}:
        return 1
    return 2


def transcript_biotype_rank(csq):
    biotype = str(csq.get("BIOTYPE", "")).lower()
    if biotype == "processed_transcript":
        return 0
    if biotype in {"nonsense_mediated_decay", "non_stop_decay", "retained_intron"}:
        return 1
    if biotype.endswith("_transcript"):
        return 1
    return 2


def consequence_above_cutoff_rank(csq, order_map, cutoff):
    cutoff_rank = order_map.get(cutoff)
    if cutoff_rank is None:
        return 1
    return 0 if consequence_rank(csq, order_map) < cutoff_rank else 1


def impactful_consequence_rank(csq, order_map):
    return consequence_above_cutoff_rank(csq, order_map, "IMPACTFUL_CUTOFF")


def genic_consequence_rank(csq, order_map):
    return consequence_above_cutoff_rank(csq, order_map, "GENIC_CUTOFF")


def weak_gene_symbol(symbol):
    if not present(symbol):
        return False
    text = str(symbol).upper()
    if text.startswith(("ENSG", "LOC", "LINC", "MIR", "RN7", "RNU", "SNORD", "SNORA", "SCARNA", "Y_RNA")):
        return True
    for prefix in ("AC", "AL", "AP"):
        if text.startswith(prefix) and len(text) > len(prefix) and text[len(prefix)].isdigit():
            return True
    return text.startswith("RP") and len(text) > 2 and (text[2].isdigit() or text[2] == "-")


def strong_gene_rank(csq):
    return 0 if not weak_gene_symbol(gene_symbol(csq)) else 1


def csq_sort_key(csq, order_map, mode):
    return (
        pseudogene_rank(csq),
        impactful_consequence_rank(csq, order_map),
        protein_coding_biotype_rank(csq),
        consequence_rank(csq, order_map),
        strong_gene_rank(csq),
        mane_rank(csq),
        canonical_rank(csq),
        transcript_biotype_rank(csq),
        genic_consequence_rank(csq, order_map),
        csq.get("Feature", ""),
        csq.get("_csq_index", 0),
    )


def choose_primary_csq(csq_records, order_map, mode):
    if not csq_records:
        return None
    return sorted(csq_records, key=lambda csq: csq_sort_key(csq, order_map, mode))[0]


def variant_key(record):
    return (record.chrom, record.pos, record.ref, record.alts[0] if record.alts else "")


def gt_string(sample_data, ref, alts):
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
    return ("|" if sample_data.phased else "/").join(alleles)


def sample_alt_depth_value(sample_data):
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


def sample_alt_depth(sample_data):
    value = sample_alt_depth_value(sample_data)
    return "" if value is None else str(value)


def sample_depth(sample_data):
    dp = sample_data.get("DP")
    return "" if dp is None else str(dp)


def sample_gq(sample_data):
    gq = sample_data.get("GQ")
    return "" if gq is None else str(gq)


def sample_format_value(sample_data, field):
    value = sample_data.get(field)
    if value is None:
        return ""
    if isinstance(value, tuple):
        return ",".join(str(item) for item in value if present(item))
    text = str(value)
    return "" if text in MISSING else text


def zygosity(sample_data, ref, alts, chrom=""):
    genotype = gt_string(sample_data, ref, alts).replace("|", "/").replace("./.", "Missing")
    if "Missing" in genotype:
        return "Missing"
    if sample_alt_depth_value(sample_data) == 0:
        return "-"
    alleles = genotype.split("/")
    if len(alleles) == 2:
        if alleles[0] == alleles[1]:
            return "-" if alleles[0] == ref else "Hom"
        return "Het"
    chrom_text = str(chrom)
    if "X" in chrom_text or "Y" in chrom_text:
        if genotype == ".":
            return "Missing"
        return "-" if genotype == ref else "Hom"
    return genotype
