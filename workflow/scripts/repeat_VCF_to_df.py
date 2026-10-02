#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

from pysam import VariantFile


# Keep the flattened ExpansionHunter output close to the columns previously emitted from DRAGEN VCFs.
OUTPUT_COLUMNS = [
    "SAMPLE",
    "STRCHIVE_LOCUS_ID",
    "CHROM",
    "POS",
    "END",
    "VARID",
    "REF",
    "RL",
    "RU",
    "GT",
    "SO",
    "REPCN",
    "REPCI",
    "ADSP",
    "ADFL",
    "ADIR",
    "LC",
    "FILTER",
]

def normalize_chromosome(chromosome):
    return chromosome if chromosome.startswith("chr") else f"chr{chromosome}"


def load_targets(threshold_path):
    # ExpansionHunter emits these stable component IDs from the same catalog.
    targets = {}
    with open(threshold_path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            targets[row["Target VariantId"]] = row["STRchive LocusId"]
    return targets


def format_value(value, precision=None):
    # Preserve multi-valued VCF fields with slash separators and represent missing values as dots.
    if value is None:
        return "."
    if isinstance(value, (tuple, list)):
        return "/".join(format_value(item, precision) for item in value)
    if isinstance(value, float) and precision is not None:
        return str(round(value, precision))
    return str(value)


def genotype(sample_call):
    alleles = sample_call.get("GT")
    if alleles is None:
        return "."
    separator = "|" if sample_call.phased else "/"
    return separator.join("." if allele is None else str(allele) for allele in alleles)


def record_variant_id(record):
    variant_id = record.info.get("VARID")
    if isinstance(variant_id, (tuple, list)):
        variant_id = variant_id[0]
    return None if variant_id is None else str(variant_id)


def call_row(sample, locus_id, record, sample_call):
    # Flatten the selected VCF INFO and FORMAT fields into one report-compatible row.
    sample_value = lambda field: format_value(sample_call.get(field))
    filters = ";".join(record.filter.keys()) or "."
    return {
        "SAMPLE": sample,
        "STRCHIVE_LOCUS_ID": locus_id,
        "CHROM": normalize_chromosome(record.chrom),
        "POS": record.pos,
        "END": record.stop,
        "VARID": format_value(record.info.get("VARID")),
        "REF": format_value(record.info.get("REF")),
        "RL": format_value(record.info.get("RL")),
        "RU": format_value(record.info.get("RU")),
        "GT": genotype(sample_call),
        "SO": sample_value("SO"),
        "REPCN": sample_value("REPCN"),
        "REPCI": sample_value("REPCI"),
        "ADSP": sample_value("ADSP"),
        "ADFL": sample_value("ADFL"),
        "ADIR": sample_value("ADIR"),
        "LC": format_value(sample_call.get("LC"), precision=2),
        "FILTER": filters,
    }


def sample_calls(vcf_path, sample, targets):
    rows = []
    with VariantFile(vcf_path) as variants:
        vcf_sample = (
            sample
            if sample in variants.header.samples
            else next(iter(variants.header.samples))
        )
        for record in variants:
            # Retain only the component configured as reportable in the threshold resource.
            locus_id = targets.get(record_variant_id(record))
            if locus_id:
                rows.append(call_row(sample, locus_id, record, record.samples[vcf_sample]))
    return rows


def main(samples_tsv, expansionhunter_dir, disease_thresholds, output_file):
    with open(samples_tsv, newline="") as handle:
        samples = [row["sample"] for row in csv.DictReader(handle, delimiter="\t")]
    targets = load_targets(disease_thresholds)
    rows = [
        row
        for sample in samples
        for row in sample_calls(Path(expansionhunter_dir) / f"{sample}.vcf", sample, targets)
    ]

    with open(output_file, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} target-component calls")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples_tsv", required=True)
    parser.add_argument("--expansionhunter_dir", required=True)
    parser.add_argument("--disease_thresholds", required=True)
    parser.add_argument("--output_file", required=True)
    args = parser.parse_args()
    main(
        args.samples_tsv,
        args.expansionhunter_dir,
        args.disease_thresholds,
        args.output_file,
    )