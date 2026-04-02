import argparse
import glob
import os
import pandas as pd
import numpy as np
from pysam import VariantFile

def recode_gt(gt):
    alleles = []
    for allele in gt:
        if allele is None:
            alleles.append(".")
        else:
            alleles.append(str(allele))
    return "/".join(alleles)

def main(repeat_vcfs, output_file):
    samples_df = pd.read_csv(samples_tsv, sep="\t", dtype=str)
    results_paths = samples_df["DRAGEN_results_dir"].tolist()
    try:
        vcf_files = [glob.glob(f"{dir}/output/*repeats.vcf.gz")[0] for dir in results_paths]
    except IndexError:
        raise ValueError(f"No repeats VCF found in {results_paths}")
    with open(output_file, "w") as f:
        for repeat_vcf in vcf_files:
            variants = VariantFile(repeat_vcf)
            samples = variants.header.samples
            for rec in variants.fetch():
                # INFO fields
                VARID = rec.info["VARID"]
                REF = rec.info["REF"]  # Number of repeat units spanned by the repeat in the reference
                RL = rec.info["RL"]  # Reference length in bp
                RU = rec.info["RU"]  # Repeat unit in the reference orientation
                # Sample fields:
                for sample in samples:
                    GT = rec.samples[sample]["GT"]  # Genotype of the sample at the variant
                    GT = recode_gt(GT)
                    SO = rec.samples[sample]["SO"]  # Type of reads that support the allele. Values can be SPANNING, FLANKING, or INREPEAT. These values indicate if the reads span, flank, or are fully contained in the repeat.
                    REPCN = rec.samples[sample]["REPCN"]  # Number of repeat units spanned by the allele
                    REPCI = rec.samples[sample]["REPCI"]  # Confidence interval for REPCN
                    ADSP = rec.samples[sample]["ADSP"]  # Number of spanning reads consistent with the allele
                    ADFL = rec.samples[sample]["ADFL"]  # Number of flanking reads consistent with the allele
                    ADIR = rec.samples[sample]["ADIR"]  # Number of in-repeat reads consistent with the allele
                    LC = round(rec.samples[sample]["LC"], 2)  # Locus Coverage
                    f.write(
                        f"{sample}\t{rec.chrom}\t{rec.pos}\t{VARID}\t{REF}\t{RL}\t{RU}\t{GT}\t{SO}\t{REPCN}\t{REPCI}\t{ADSP}\t{ADFL}\t{ADIR}\t{LC}\n"
                    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Converts DRAGEN repeat VCF(s) to a TSV"
    )
    parser.add_argument(
        "--samples_tsv", type=str, help="Samples TSV with DRAGEN results paths", required=True
    )
    parser.add_argument(
        "--output_file", type=str, help="Output filename", required=True
    )

    args = parser.parse_args()
    samples_tsv = args.samples_tsv
    output_file = args.output_file

    main(samples_tsv, output_file)