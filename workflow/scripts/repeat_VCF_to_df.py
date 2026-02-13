import argparse
import glob
import os
import pandas as pd
import numpy as np
from pysam import VariantFile


def main(repeat_vcf, output_file):
    with open(output_file, "w") as f:
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
                GT = rec.samples[sample]["GT"]
                SO = rec.samples[sample]["SO"]  # Type of reads that support the allele. Values can be SPANNING, FLANKING, or INREPEAT. These values indicate if the reads span, flank, or are fully contained in the repeat.
                REPCN = rec.samples[sample]["REPCN"]  # Number of repeat units spanned by the allele
                REPCI = rec.samples[sample]["REPCI"]  # Confidence interval for REPCN
                ADSP = rec.samples[sample]["ADSP"]  # Number of spanning reads consistent with the allele
                ADFL = rec.samples[sample]["ADFL"]  # Number of flanking reads consistent with the allele
                ADIR = rec.samples[sample]["ADIR"]  # Number of in-repeat reads consistent with the allele
                LC = rec.samples[sample]["LC"]  # Locus Coverage
                f.write(
                    f"{sample}\t{rec.chrom}\t{rec.pos}\t{VARID}\t{REF}\t{RL}\t{RU}\t{GT}\t{SO}\t{REPCN}\t{REPCI}\t{ADSP}\t{ADFL}\t{ADIR}\t{LC}\n"
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Converts DRAGEN repeat VCF(s) to a TSV"
    )
    parser.add_argument(
        "--repeat_vcf", type=str, help="DRAGEN repeat VCF", required=True
    )
    parser.add_argument(
        "--output_file", type=str, help="Output filename", required=True
    )

    args = parser.parse_args()
    repeat_vcf = args.repeat_vcf
    output_file = args.output_file

    main(repeat_vcf, output_file)
