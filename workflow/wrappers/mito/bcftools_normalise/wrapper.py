#Please note that the AF field has been duplicated and named VAF; 
#Any "." or missing value in the original AF field will be changed to 0.0 in the new VAF field.
#The re-formatting of the VAF field is necessary for the downstream mity report.py script but the original AF field remains unchanged. 

from snakemake.shell import shell


log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

family = snakemake.wildcards.family
outdir = snakemake.params.outdir
tool = snakemake.params.tool
cphi_dragen_anno = snakemake.params.cphi_dragen_anno
input_vcf = snakemake.input[0]
pythonpath = tool.replace("bin", "")
python = " export PYTHONPATH={pythonpath}; "

# use mity's hg38 reference fasta 
reference_fasta = pythonpath + "/mitylib/reference/hg38.chrM.fa"

bcf_normalise = " bcftools norm -f {reference_fasta} -m-both -Oz -o {outdir}/{family}.normalise.vcf.gz {input_vcf};"
vt = " vt decompose_blocksub -o {outdir}/{family}.normalise.decompose.unformatted.vcf.gz {outdir}/{family}.normalise.vcf.gz;"
add_VAF_field = " bcftools +fill-tags {outdir}/{family}.normalise.decompose.unformatted.vcf.gz -Ou -- -t FORMAT/VAF | bcftools view -Oz -o {outdir}/{family}.normalise.decompose.temp.vcf.gz;"
reformat_empty_VAF = " python {cphi_dragen_anno}/workflow/scripts/format_missing_vaf.py {outdir}/{family}.normalise.decompose.temp.vcf.gz {outdir}/{family}.mt.normalise.decompose.vcf.gz;"
tabix = " tabix {outdir}/{family}.mt.normalise.decompose.vcf.gz;"
remove_intermediate_files = " rm {outdir}/{family}.normalise.decompose.unformatted.vcf.gz {outdir}/{family}.normalise.decompose.temp.vcf.gz {outdir}/{family}.normalise.vcf.gz"


shell("(" + python + bcf_normalise + vt + add_VAF_field + reformat_empty_VAF + tabix + remove_intermediate_files +") {log}")