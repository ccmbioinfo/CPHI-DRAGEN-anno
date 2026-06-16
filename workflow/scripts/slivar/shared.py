#!/usr/bin/env python3

import csv
from collections import defaultdict


MISSING = {"", ".", "None", "NA"}


def present(value):
    if value is None:
        return False
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
    for field in ("SYMBOL", "HGNC"):
        value = csq.get(field, "")
        if present(value):
            return value
    gene = csq.get("Gene", "")
    if present(gene) and not str(gene).startswith("ENSG"):
        return gene
    return ""


def variant_key(record):
    return (record.chrom, record.pos, record.ref, record.alts[0] if record.alts else "")


def gt_string(sample_data, ref, alts):
    gt = sample_data.get("GT")
    if gt is None:
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


def sample_alt_depth(sample_data):
    ad = sample_data.get("AD")
    if ad is None or len(ad) < 2:
        return ""
    try:
        return str(ad[1])
    except Exception:
        return ""


def sample_depth(sample_data):
    dp = sample_data.get("DP")
    return "" if dp is None else str(dp)


def sample_gq(sample_data):
    gq = sample_data.get("GQ")
    return "" if gq is None else str(gq)


def zygosity(sample_data, chrom=""):
    gt = sample_data.get("GT")
    ad = sample_data.get("AD")
    if gt is None or all(allele is None for allele in gt):
        return "Missing"
    alt_ad = 0
    if ad is not None and len(ad) > 1 and ad[1] is not None:
        try:
            alt_ad = int(ad[1])
        except Exception:
            alt_ad = 0
    if alt_ad == 0:
        return "-"
    called = [allele for allele in gt if allele is not None]
    if not called:
        return "Missing"
    if all(allele == 0 for allele in called):
        return "-"
    if len(called) == 1 and called[0] != 0:
        chrom_text = str(chrom)
        if "X" in chrom_text or "Y" in chrom_text:
            return "Hom"
    if len(called) == 2 and called[0] == called[1] and called[0] != 0:
        return "Hom"
    return "Het"
