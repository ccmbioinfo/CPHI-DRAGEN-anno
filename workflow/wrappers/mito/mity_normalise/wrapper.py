#this wrapper is not used in the pipeline, the rule has been replaced by bcftools_normalise 

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True,append=True)

family = snakemake.wildcards.family
outdir = snakemake.params.outdir
tool = snakemake.params.tool
pythonpath = tool.replace("bin", "")
reference = "hg38"

python = " export PYTHONPATH={pythonpath}; "
mity = " {tool}/mity normalise --reference {reference} --prefix {family}.mity --output-dir {outdir} {snakemake.input};"
vt = " vt decompose_blocksub -o {outdir}/{family}.mity.normalise.decompose.vcf.gz {outdir}/{family}.mity.normalise.vcf.gz;"
tabix = " tabix {outdir}/{family}.mity.normalise.decompose.vcf.gz "
shell("(" + python + mity + vt + tabix + ") {log}")