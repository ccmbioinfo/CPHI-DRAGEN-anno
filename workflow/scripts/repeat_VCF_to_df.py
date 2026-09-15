import argparse
import json
import re

import pandas as pd
from pysam import VariantFile


OUTPUT_COLUMNS = ["SAMPLE","CHROM","POS","VARID","REF","RL","RU","GT","SO","REPCN","REPCI","ADSP","ADFL","ADIR","LC",]

# Native DRAGEN VARIDs that differ from the threshold report keys.
DRAGEN_GENE_NAMES = {
    "HOXA13_1": "HOXA13-I",
    "HOXA13_2": "HOXA13-II",
    "HOXA13_3": "HOXA13-III",
    "ARX_1": "EIEE1_ARX",
    "ARX_2": "PRTS_ARX",
    "C9ORF72": "C9orf72",
}

# Catalog IDs whose report key cannot be taken from the final underscore field.
CATALOG_REPORT_NAMES = {
    "pre-MIR7-2_CHNG3": "pre-MIR7-2",
}


def recode_gt(gt):
    alleles = []
    for allele in gt:
        if allele is None:
            alleles.append(".")
        else:
            alleles.append(str(allele))
    return "/".join(alleles)


def vcf_to_df(vcf_file):
    rows = []
    with VariantFile(vcf_file) as variants:
        for rec in variants:
            # INFO fields
            VARID = rec.info["VARID"]
            REF = rec.info["REF"]  # Number of repeat units in the reference
            RL = rec.info["RL"]  # Reference length in bp
            RU = rec.info["RU"]  # Repeat unit in the reference orientation

            # Sample fields
            for sample in variants.header.samples:
                GT = recode_gt(rec.samples[sample]["GT"])  # Genotype
                SO = rec.samples[sample]["SO"]  # Supporting-read types
                REPCN = rec.samples[sample]["REPCN"]  # Allele repeat counts
                REPCI = rec.samples[sample]["REPCI"]  # REPCN confidence intervals
                ADSP = rec.samples[sample]["ADSP"]  # Spanning-read support
                ADFL = rec.samples[sample]["ADFL"]  # Flanking-read support
                ADIR = rec.samples[sample]["ADIR"]  # In-repeat-read support
                LC = round(rec.samples[sample]["LC"], 2)  # Locus coverage
                rows.append(
                    [
                        sample,
                        rec.chrom,
                        rec.pos,
                        VARID,
                        REF,
                        RL,
                        RU,
                        GT,
                        SO,
                        REPCN,
                        REPCI,
                        ADSP,
                        ADFL,
                        ADIR,
                        LC,
                        rec.stop,
                    ]
                )
    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS + ["END"])


def parse_region(region):
    chrom, coordinates = region.split(":", 1)
    start, end = coordinates.split("-", 1)
    return chrom, int(start), int(end)


def get_catalog_targets(variant_catalog, disease_thresholds):
    # The threshold table supplies the report key and identifies which repeat
    # component should be reported for each disease locus.
    thresholds = pd.read_csv(disease_thresholds, sep="\t", dtype=str)
    thresholds = thresholds.set_index("Gene", drop=False)
    with open(variant_catalog) as handle:
        catalog = json.load(handle)

    targets = []
    for locus in catalog:
        locus_id = locus["LocusId"]
        gene = CATALOG_REPORT_NAMES.get(locus_id)
        if gene is None:
            # Keep IDs that are already threshold keys; otherwise extract the
            # gene from STRchive's disease_gene locus naming convention.
            gene = (
                locus_id
                if locus_id in thresholds.index
                else locus_id.rsplit("_", 1)[-1]
            )
        threshold = thresholds.loc[gene]

        # ExpansionHunter represents compound loci with parallel lists of
        # reference regions and repeat motifs; simple loci contain one string.
        regions = locus["ReferenceRegion"]
        if isinstance(regions, str):
            regions = [regions]
        motifs = re.findall(r"\(([ACGTN]+)\)\*", locus["LocusStructure"], re.I)

        # Threshold coordinates are 1-based, while catalog region starts and
        # ExpansionHunter VCF positions use the catalog's 0-based start value.
        threshold_region = None
        if pd.notna(threshold["Coordinates (hg38)"]):
            threshold_region = parse_region(threshold["Coordinates (hg38)"])
        threshold_motif = None
        if pd.notna(threshold["Motif (reference orientation)"]):
            threshold_motif = threshold["Motif (reference orientation)"].upper()

        components = []
        for region, motif in zip(regions, motifs):
            chrom, pos, end = parse_region(region)
            components.append(
                {
                    "CHROM": chrom,
                    "POS": pos,
                    "END": end,
                    "GENE": gene,
                    "RU": motif.upper(),
                }
            )

        # Prefer the component identified by the threshold coordinates. Some
        # thresholds span a whole compound locus, so use its motif next. The
        # final fallback handles single-component legacy loci.
        target = next(
            (
                component
                for component in components
                if (component["CHROM"], component["POS"] + 1, component["END"])
                == threshold_region
            ),
            None,
        )
        if target is None and threshold_motif:
            target = next(
                (
                    component
                    for component in components
                    if component["RU"] == threshold_motif
                ),
                None,
            )
        if target is None:
            target = components[0]
        # Do not make a disease prediction for compound loci because
        # their motif counts cannot be interpreted by the simple threshold.
        target["MULTI_MOTIF"] = len(components) > 1
        targets.append(target)

    return pd.DataFrame(targets)


def exact_catalog_matches(repeats, catalog_targets):
    matches = repeats.merge(catalog_targets, on=["CHROM", "POS", "END", "RU"], how="inner",)
    matches["VARID"] = matches["GENE"]
    return matches[OUTPUT_COLUMNS + ["END", "GENE"]]


def remove_covered(repeats, covered):
    sample_genes = zip(repeats["SAMPLE"], repeats["GENE"])
    return repeats[[pair not in covered for pair in sample_genes]].copy()


def main(samples_tsv, expansionhunter_dir, variant_catalog, disease_thresholds, output_file,):
    samples = pd.read_csv(samples_tsv, sep="\t", dtype=str)
    catalog_targets = get_catalog_targets(variant_catalog, disease_thresholds)

    dragen_repeats = pd.concat([vcf_to_df(vcf) for vcf in samples["STR"]], ignore_index=True)
    expansionhunter_repeats = pd.concat([vcf_to_df(f"{expansionhunter_dir}/{sample}.vcf") for sample in samples["sample"]], ignore_index=True,)

    # Keep the disease-labelled DRAGEN calls, as in the original workflow.
    is_unlabelled = dragen_repeats["VARID"].str.startswith("chr")
    final_repeats = dragen_repeats[~is_unlabelled].copy()
    dragen_to_match = dragen_repeats[is_unlabelled].copy()
    final_repeats["GENE"] = final_repeats["VARID"].replace(DRAGEN_GENE_NAMES)

    # Add exact catalog matches from the remaining unlabelled DRAGEN calls.
    matched_dragen = exact_catalog_matches(dragen_to_match, catalog_targets)
    covered = set(zip(final_repeats["SAMPLE"], final_repeats["GENE"]))
    matched_dragen = remove_covered(matched_dragen, covered)
    final_repeats = pd.concat([final_repeats, matched_dragen], ignore_index=True)

    # Keep one reportable component per EH locus, then remove loci already
    # supplied by DRAGEN and append what remains.
    expansionhunter_repeats = exact_catalog_matches(expansionhunter_repeats, catalog_targets)
    covered = set(zip(final_repeats["SAMPLE"], final_repeats["GENE"]))
    expansionhunter_repeats = remove_covered(expansionhunter_repeats, covered)
    final_repeats = pd.concat([final_repeats, expansionhunter_repeats], ignore_index=True)

    multi_motif_genes = set(
        catalog_targets.loc[catalog_targets["MULTI_MOTIF"], "GENE"]
    )
    final_repeats["MULTI_MOTIF"] = final_repeats["GENE"].isin(multi_motif_genes)
    final_repeats[OUTPUT_COLUMNS + ["MULTI_MOTIF"]].to_csv(
        output_file, sep="\t", index=False, header=False
    )

    print(f"DRAGEN disease-labelled calls: {sum(~is_unlabelled)}")
    print(f"DRAGEN exact catalog matches added: {len(matched_dragen)}")
    print(f"ExpansionHunter fallback calls added: {len(expansionhunter_repeats)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combines DRAGEN and ExpansionHunter repeat calls")
    parser.add_argument("--samples_tsv", required=True)
    parser.add_argument("--family", required=True)
    parser.add_argument("--expansionhunter_dir", required=True)
    parser.add_argument("--variant_catalog", required=True)
    parser.add_argument("--disease_thresholds", required=True)
    parser.add_argument("--output_file", required=True)
    args = parser.parse_args()

    main(
        args.samples_tsv,
        args.expansionhunter_dir,
        args.variant_catalog,
        args.disease_thresholds,
        args.output_file,
    )