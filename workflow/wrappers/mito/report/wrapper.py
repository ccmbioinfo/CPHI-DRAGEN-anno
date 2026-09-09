from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True,append=True)

family = snakemake.wildcards.family
outdir = snakemake.params.outdir

chr_contig = "chrM"
base_path= snakemake.params.base_path
vcfanno_config= snakemake.params.vcfanno_config
report_config= snakemake.params.report_config

mity = " mity report -k --contig {chr_contig} --prefix {family} --vcfanno-base-path {base_path} --custom-vcfanno-config {vcfanno_config} --custom-report-config {report_config} --output-dir {outdir} {snakemake.input}; "
rename_vcf_file= " mv mitochondrial_variants/{family}.mt.normalise.decompose.mity.annotated.vcf mitochondrial_variants/{family}.mity.annotated.vcf ; "
compress_vcf = " bgzip mitochondrial_variants/{family}.mity.annotated.vcf"

shell("(" + mity + rename_vcf_file + compress_vcf + ") {log}")