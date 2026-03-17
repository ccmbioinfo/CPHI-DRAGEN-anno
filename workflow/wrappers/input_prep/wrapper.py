from snakemake.shell import shell
import pandas as pd

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
family = snakemake.wildcards.family
units = pd.read_table(snakemake.input.units, dtype=str).set_index(["family"], drop=False)
input_vcf=units.loc[family, "small_variant_vcf"]
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "( bcftools filter -e 'STRLEN(REF) >=50 | STRLEN(ALT) >=50' {input_vcf} -O z -o filtered/{family}.vcf.gz ) {log}"
)