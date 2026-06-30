#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def load_keys(path):
    keys = {}

    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)

        for row in reader:
            position = (row.get("Position") or "").strip()
            ref = (row.get("Ref") or "").strip()
            alt = (row.get("Alt") or "").strip()

            if not position or not ref or not alt:
                continue

            key = f"{position}:{ref}:{alt}"
            keys[key] = row

    return keys


def write_key_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["key"])
        for key in rows:
            writer.writerow([key])


def main():
    parser = argparse.ArgumentParser(
        description="Compare Gemini and slivar final report CSVs by Position/Ref/Alt key."
    )
    parser.add_argument("--gemini-report", required=True)
    parser.add_argument("--slivar-report", required=True)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    gemini_keys = load_keys(args.gemini_report)
    slivar_keys = load_keys(args.slivar_report)

    gemini_set = set(gemini_keys)
    slivar_set = set(slivar_keys)

    shared = sorted(gemini_set & slivar_set)
    gemini_only = sorted(gemini_set - slivar_set)
    slivar_only = sorted(slivar_set - gemini_set)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    summary_path = Path(f"{out_prefix}.summary.tsv")
    with open(summary_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "count"])
        writer.writerow(["gemini_total", len(gemini_set)])
        writer.writerow(["slivar_total", len(slivar_set)])
        writer.writerow(["shared", len(shared)])
        writer.writerow(["gemini_only", len(gemini_only)])
        writer.writerow(["slivar_only", len(slivar_only)])

    write_key_csv(Path(f"{out_prefix}.shared.csv"), shared)
    write_key_csv(Path(f"{out_prefix}.gemini_only.csv"), gemini_only)
    write_key_csv(Path(f"{out_prefix}.slivar_only.csv"), slivar_only)


if __name__ == "__main__":
    main()
