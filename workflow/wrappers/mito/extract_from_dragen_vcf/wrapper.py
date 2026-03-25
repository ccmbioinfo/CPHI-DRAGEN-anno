from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

family = snakemake.wildcards.family

dragen_vcf = snakemake.input.vcf

out = snakemake.output[0]


# Extract mitochondrial variants from the Dragen VCF and output as a new VCF

shell(
    """
    bcftools index -f -t {dragen_vcf}
    bcftools view -r chrM {dragen_vcf} -Oz -o {out}
    bcftools index -f -t {out}
    {log}
    """
)