#!/usr/bin/env python3

import argparse
import csv
from datetime import date
import os


CATEGORIES = [
    ("BENIGN", "Benign range(s)"),
    ("INTERMEDIATE", "Intermediate range(s)"),
    ("PATHOGENIC", "Pathogenic range(s)"),
]
CALL_FIELDS = [
    ("GT", "GT"),
    ("motif_count", "REPCN"),
    ("REPCI", "REPCI"),
    ("FILTER", "FILTER"),
    ("SO", "SO"),
    ("ADSP", "ADSP"),
    ("ADFL", "ADFL"),
    ("ADIR", "ADIR"),
    ("coverage", "LC"),
]


def ranges(value):
    parsed = []
    for interval in filter(None, value.split(";")):
        minimum, separator, maximum = interval.partition("-")
        parsed.append(
            (int(minimum), None if separator and maximum == "*" else int(maximum or minimum))
        )
    return parsed


def classify_allele(count, threshold):
    for category, field in CATEGORIES:
        for minimum, maximum in ranges(threshold[field]):
            if count >= minimum and (maximum is None or count <= maximum):
                return category
    return "UNKNOWN"


def classify(repcn, threshold):
    try:
        counts = [int(float(value)) for value in repcn.replace("|", "/").split("/")]
    except (AttributeError, ValueError):
        return "MISSING"
    if threshold["Classification enabled"].lower() != "true":
        return "UNKNOWN"

    allele_classes = [classify_allele(count, threshold) for count in counts]
    if "PATHOGENIC" in allele_classes:
        return "PATHOGENIC"
    if "UNKNOWN" in allele_classes:
        return "UNKNOWN"
    if "INTERMEDIATE" in allele_classes:
        return "INTERMEDIATE"
    return "BENIGN"


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def first_call(calls, samples, field):
    for sample in samples:
        value = calls.get(sample, {}).get(field, ".")
        if value not in ("", "."):
            return value
    return "."


def build_report(repeat_tsv, thresholds_tsv, samples_tsv):
    samples = [row["sample"] for row in read_tsv(samples_tsv)]
    calls = {
        (row["STRCHIVE_LOCUS_ID"], row["SAMPLE"]): row
        for row in read_tsv(repeat_tsv)
    }
    report = []

    for threshold in read_tsv(thresholds_tsv):
        locus_id = threshold["STRchive LocusId"]
        locus_calls = {sample: calls.get((locus_id, sample), {}) for sample in samples}
        chromosome, coordinates = threshold["Target region (hg38)"].split(":", 1)
        start = coordinates.split("-", 1)[0]
        row = {
            "CHROM": chromosome,
            "POS": start,
            "REF_REPEAT_COUNT": first_call(locus_calls, samples, "REF"),
            "REF_LEN_BP": first_call(locus_calls, samples, "RL"),
            "MOTIF": threshold["Target motif"],
            "GENE": threshold["Gene"],
            "DISORDER": threshold["Disorder"],
            "DISEASE_THRESHOLD": threshold["Disease threshold"] or ".",
        }
        for sample in samples:
            row[f"{sample}_DISEASE_PREDICTION"] = classify(
                locus_calls[sample].get("REPCN", "."), threshold
            )
        for sample in samples:
            for output_name, input_name in CALL_FIELDS:
                row[f"{sample}_{output_name}"] = (
                    locus_calls[sample].get(input_name, ".") or "."
                )
        row.update(
            {
                "STRCHIVE_LOCUS_ID": locus_id,
                "BENIGN_RANGES": threshold["Benign range(s)"] or ".",
                "INTERMEDIATE_RANGES": threshold["Intermediate range(s)"] or ".",
                "PATHOGENIC_RANGES": threshold["Pathogenic range(s)"] or ".",
                "TARGET_REGION": threshold["Target region (hg38)"],
                "TARGET_VARIANT_ID": threshold["Target VariantId"],
                "STRCHIVE_URL": threshold["STRchive URL"],
            }
        )
        report.append(row)
    return report, samples


def main(repeat_tsv, disease_thresholds, samples_tsv, output_file):
    rows, samples = build_report(repeat_tsv, disease_thresholds, samples_tsv)
    leading = [
        "CHROM",
        "POS",
        "REF_REPEAT_COUNT",
        "REF_LEN_BP",
        "MOTIF",
        "GENE",
        "DISORDER",
        "DISEASE_THRESHOLD",
    ]
    predictions = [f"{sample}_DISEASE_PREDICTION" for sample in samples]
    sample_fields = [
        f"{sample}_{output_name}"
        for sample in samples
        for output_name, _ in CALL_FIELDS
    ]
    resource_fields = [
        "STRCHIVE_LOCUS_ID",
        "BENIGN_RANGES",
        "INTERMEDIATE_RANGES",
        "PATHOGENIC_RANGES",
        "TARGET_REGION",
        "TARGET_VARIANT_ID",
        "STRCHIVE_URL",
    ]

    output_prefix = output_file.removesuffix(".hg38.csv")
    dated_output = f"{output_prefix}.{date.today().isoformat()}.hg38.csv"
    os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
    with open(dated_output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=leading + predictions + sample_fields + resource_fields,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    if os.path.lexists(output_file):
        os.remove(output_file)
    os.symlink(os.path.basename(dated_output), output_file)
    print(f"Wrote {dated_output} ({len(rows)} loci)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat_tsv", required=True)
    parser.add_argument("--disease_thresholds", required=True)
    parser.add_argument("--samples_tsv", required=True)
    parser.add_argument("--output_file", required=True)
    args = parser.parse_args()
    main(
        args.repeat_tsv,
        args.disease_thresholds,
        args.samples_tsv,
        args.output_file,
    )
