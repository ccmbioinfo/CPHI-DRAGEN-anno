from snakemake.shell import shell

cram = snakemake.input.cram
out_prefix = snakemake.params.out_prefix
selfsm = snakemake.output.selfsm
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
sample = snakemake.params.sample
ref = snakemake.params.ref
svdp = snakemake.params.svdp

shell(
    """
    verifybamid2 \
        --BamFile {cram} \
        --Reference {ref} \
        --SVDPrefix {svdp} \
        --Output {out_prefix} \
        {log}

    awk -v sample="{sample}" '
        BEGIN {{ OFS="\t" }}
        /^#/ {{ print; next }}
        $1=="DefaultSampleName" {{ $1=sample }}
        {{ print }}
    ' {out_prefix}.selfSM > {out_prefix}.selfSM.tmp

    mv {out_prefix}.selfSM.tmp {selfsm}
    """
)
