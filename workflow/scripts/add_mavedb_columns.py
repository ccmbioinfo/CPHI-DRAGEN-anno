import argparse
import logging
import os
from collections import defaultdict
from typing import Optional

import pandas as pd


DEFAULT_REPORT_TYPES = [
    "wgs.coding",
    "wgs.high.impact",
    "wgs.denovo",
    "panel",
    "panel-flank",
]

MAVEDB_COLUMNS = [
    "MaveDB_link",
    "MaveDB_assertion",
    "MaveDB_functional_class",
    "MaveDB_score",
]


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def clean_value(value) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def is_missing(value) -> bool:
    return clean_value(value) in {"", "NA", ".", "None", "nan"}

# Preserve first-seen order so collapsed annotations stay stable across runs.
def deduplicate_preserve_order(values):
    seen = set()
    ordered = []
    for value in values:
        if value not in seen:
            seen.add(value)
            ordered.append(value)
    return ordered


def collapse_values(values: list[str]) -> str:
    filtered = [value for value in values if not is_missing(value)]
    filtered = deduplicate_preserve_order(filtered)
    if not filtered:
        return "."
    return ";".join(filtered)


def format_mavedb_link(url: str) -> str:
    if is_missing(url):
        return "."
    label = "MAVE_variant_link" if "/variants/" in url else "MAVE_experiment_link"
    return f'=HYPERLINK("{url}","{label}")'

# Choose one link per MaveDB row, preferring a variant page over broader experiment pages.
def pick_row_link(variant_url: str, preferred_link: str, scoreset_url: str) -> str:
    for candidate in (variant_url, preferred_link, scoreset_url):
        if not is_missing(candidate):
            return format_mavedb_link(clean_value(candidate))
    return "."

# If multiple MaveDB rows match, keep a single final link and prefer a variant-level page.
def collapse_link(values: list[str]) -> str:
    filtered = [value for value in values if not is_missing(value)]
    filtered = deduplicate_preserve_order(filtered)
    if not filtered:
        return "."

    variant_links = [value for value in filtered if "MAVE_variant_link" in value]
    if variant_links:
        return variant_links[0]

    experiment_links = [value for value in filtered if "MAVE_experiment_link" in value]
    if experiment_links:
        return experiment_links[0]

    return filtered[0]

# Parse the CRE Info field and build NM:p. and NP:p. fallback keys from any transcript/protein entries present.
def extract_info_match_keys(info_value: str) -> tuple[list[str], list[str]]:
    nm_keys = []
    np_keys = []

    if is_missing(info_value):
        return nm_keys, np_keys

    for unit in clean_value(info_value).split(","):
        if is_missing(unit):
            continue

        parts = [clean_value(part) for part in unit.split(":") if clean_value(part)]
        protein_change = next((part for part in parts if part.startswith("p.")), "")
        if not protein_change:
            continue

        nm_ids = [part for part in parts if part.startswith("NM_")]
        np_ids = [part for part in parts if part.startswith("NP_")]

        nm_keys.extend([f"{nm_id}:{protein_change}" for nm_id in nm_ids])
        np_keys.extend([f"{np_id}:{protein_change}" for np_id in np_ids])

    return deduplicate_preserve_order(nm_keys), deduplicate_preserve_order(np_keys)


def extract_refseq_key(refseq_change: str) -> str:
    if is_missing(refseq_change):
        return ""

    parts = [clean_value(part) for part in clean_value(refseq_change).split(":") if clean_value(part)]
    if not parts:
        return ""

    accession = next((part for part in parts if part.startswith("NM_")), "")
    protein_change = next((part for part in parts if part.startswith("p.")), "")
    if accession and protein_change:
        return f"{accession}:{protein_change}"
    return ""

# Collapse all matched MaveDB rows for one lookup key into the four reportable output fields.
def collect_match(records: list[dict]) -> dict[str, str]:
    if not records:
        return {column: "." for column in MAVEDB_COLUMNS}

    return {
        "MaveDB_link": collapse_link([record["MaveDB_link"] for record in records]),
        "MaveDB_assertion": collapse_values([record["MaveDB_assertion"] for record in records]),
        "MaveDB_functional_class": collapse_values([record["MaveDB_functional_class"] for record in records]),
        "MaveDB_score": collapse_values([record["MaveDB_score"] for record in records]),
    }


def has_annotation(match: dict[str, str]) -> bool:
    return any(not is_missing(match[column]) for column in MAVEDB_COLUMNS)

# Return the first lookup key that yields any usable MaveDB annotation.
def first_match_for_keys(keys: list[str], lookup: dict[str, list[dict]]) -> Optional[dict[str, str]]:
    for key in keys:
        records = lookup.get(key, [])
        if records:
            match = collect_match(records)
            if has_annotation(match):
                return match
    return None

# Load the compact MaveDB resource once and build separate lookup tables for genomic, NM-based, and NP-based matching.
def load_mavedb(mavedb_tsv: str) -> tuple[dict[str, list[dict]], dict[str, list[dict]], dict[str, list[dict]]]:
    if not os.path.exists(mavedb_tsv):
        raise FileNotFoundError(f"MaveDB file not found: {mavedb_tsv}")

    header = pd.read_csv(mavedb_tsv, sep="\t", nrows=0)
    available_columns = set(header.columns)
    # Only load the fields needed for report annotation and fallback matching.
    desired_columns = [
        "mavedb_preferred_link",
        "mavedb_variant_url",
        "mavedb_scoreset_url",
        "display_assertion",
        "threshold_functional_class",
        "score",
        "nm_accessions",
        "np_accession",
        "hgvs_pro",
        "Position",
        "Ref",
        "Alt",
    ]
    usecols = [column for column in desired_columns if column in available_columns]

    mavedb_df = pd.read_csv(
        mavedb_tsv,
        sep="\t",
        usecols=usecols,
        dtype=str,
        keep_default_na=False,
    )
    
    # Backfill absent columns so the script can tolerate slightly different resource versions.
    for column in desired_columns:
        if column not in mavedb_df.columns:
            mavedb_df[column] = ""

    mavedb_df["MaveDB_link"] = mavedb_df.apply(
        lambda row: pick_row_link(
            row["mavedb_variant_url"],
            row["mavedb_preferred_link"],
            row["mavedb_scoreset_url"],
        ),
        axis=1,
    )
    mavedb_df["MaveDB_assertion"] = mavedb_df["display_assertion"].map(
        lambda value: clean_value(value) if not is_missing(value) else "."
    )
    mavedb_df["MaveDB_functional_class"] = mavedb_df["threshold_functional_class"].map(
        lambda value: clean_value(value) if not is_missing(value) else "."
    )
    mavedb_df["MaveDB_score"] = mavedb_df["score"].map(
        lambda value: clean_value(value) if not is_missing(value) else "."
    )

    # Build separate indexes for genomic, RefSeq transcript, and RefSeq protein matching.
    genomic_lookup = defaultdict(list)
    nm_lookup = defaultdict(list)
    np_lookup = defaultdict(list)

    for row in mavedb_df.to_dict(orient="records"):
        record = {
            "MaveDB_link": row["MaveDB_link"],
            "MaveDB_assertion": row["MaveDB_assertion"],
            "MaveDB_functional_class": row["MaveDB_functional_class"],
            "MaveDB_score": row["MaveDB_score"],
        }

        position = clean_value(row["Position"])
        ref = clean_value(row["Ref"])
        alt = clean_value(row["Alt"])
        if position and ref and alt:
            genomic_lookup[f"{position}-{ref}-{alt}"].append(record)

        nm_accession = clean_value(row["nm_accessions"])
        hgvs_pro = clean_value(row["hgvs_pro"])
        if nm_accession and hgvs_pro:
            nm_lookup[f"{nm_accession}:{hgvs_pro}"].append(record)

        np_accession = clean_value(row["np_accession"])
        if np_accession and hgvs_pro:
            np_lookup[f"{np_accession}:{hgvs_pro}"].append(record)

    logging.info("Loaded %d MaveDB rows from %s", len(mavedb_df), mavedb_tsv)
    return genomic_lookup, nm_lookup, np_lookup


def annotate_report(
    report_path: str,
    genomic_lookup: dict[str, list[dict]],
    nm_lookup: dict[str, list[dict]],
    np_lookup: dict[str, list[dict]],
) -> None:
    report_df = pd.read_csv(report_path, dtype=str, keep_default_na=False)
    logging.info("Annotating %s (%d rows)", report_path, len(report_df))

    for column in MAVEDB_COLUMNS:
        report_df[column] = "."

    for index, row in report_df.iterrows():
        position = clean_value(row.get("Position", ""))
        ref = clean_value(row.get("Ref", ""))
        alt = clean_value(row.get("Alt", ""))
        # Match in priority order: genomic coordinates first, then RefSeq protein change, then Info-derived NM/NP protein keys.
        matched = None
        if position and ref and alt:
            matched = collect_match(genomic_lookup.get(f"{position}-{ref}-{alt}", []))
            if not has_annotation(matched):
                matched = None

        if matched is None:
            refseq_key = extract_refseq_key(row.get("Refseq_change", ""))
            if refseq_key:
                matched = collect_match(nm_lookup.get(refseq_key, []))
                if not has_annotation(matched):
                    matched = None

        if matched is None:
            info_nm_keys, info_np_keys = extract_info_match_keys(row.get("Info", ""))
            matched = first_match_for_keys(info_nm_keys, nm_lookup)
            if matched is None:
                matched = first_match_for_keys(info_np_keys, np_lookup)

        if matched is None:
            continue

        for column in MAVEDB_COLUMNS:
            report_df.at[index, column] = matched[column]
    
    # Keep the MaveDB columns near SpliceAI so they stay with the other functional impact annotations.
    existing_columns = [column for column in report_df.columns if column not in MAVEDB_COLUMNS]
    if "SpliceAI_score" in existing_columns:
        insert_at = existing_columns.index("SpliceAI_score") + 1
        ordered_columns = existing_columns[:insert_at] + MAVEDB_COLUMNS + existing_columns[insert_at:]
    else:
        ordered_columns = existing_columns + MAVEDB_COLUMNS
    report_df = report_df[ordered_columns]

    # Resolve the report symlink so we rewrite the dated backing CSV rather than the symlink itself.
    target_path = os.path.realpath(report_path) if os.path.islink(report_path) else report_path
    report_df.to_csv(target_path, index=False)
    logging.info("Wrote annotated report to %s", target_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Add MaveDB columns to report CSVs in place.")
    parser.add_argument("--family", required=True)
    parser.add_argument("--mavedb-tsv", required=True)
    parser.add_argument("--reports-dir", default="reports")
    parser.add_argument("--suffix", default="hg38.csv")
    parser.add_argument("--report-types", nargs="+", default=DEFAULT_REPORT_TYPES)
    args = parser.parse_args()

    configure_logging()

    genomic_lookup, nm_lookup, np_lookup = load_mavedb(args.mavedb_tsv)

    annotated_count = 0
    for report_type in args.report_types:
        report_path = os.path.join(args.reports_dir, f"{args.family}.{report_type}.CH.{args.suffix}")
        if not os.path.exists(report_path):
            logging.info("Skipping missing report %s", report_path)
            continue

        annotate_report(report_path, genomic_lookup, nm_lookup, np_lookup)
        annotated_count += 1

    logging.info("Annotated %d report(s) for family %s", annotated_count, args.family)


if __name__ == "__main__":
    main()
