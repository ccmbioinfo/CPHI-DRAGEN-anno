rule cram_to_bam:
    input:
        cram=get_cram
    output:
        bam=temp("qc/qualimap/{family}/{sample}.bam")
    log:
        "logs/qc/qualimap/{family}/{sample}.cram_to_bam.log"
    conda:
        "../envs/samtools.yaml"
    threads: 8
    params:
        ref=config["ref"]["genome"]
    shell:
        """
        mkdir -p qc/qualimap/{wildcards.family} logs/qc/qualimap/{wildcards.family}
        samtools view -@ {threads} -T {params.ref} -b {input.cram} > {output.bam}
        """

rule peddy:
    input:
        vcf=get_sequence_var_vcf,
        ped=config["run"]["ped"]
    output:
        pca="qc/peddy/{family}.background_pca.json",
        pca_png="qc/peddy/{family}.pca_check.png",
        html="qc/peddy/{family}.html",
        het="qc/peddy/{family}.het_check.csv",
        het_png="qc/peddy/{family}.het_check.png",
        vs="qc/peddy/{family}.vs.html",
        sex="qc/peddy/{family}.sex_check.csv",
        sex_png="qc/peddy/{family}.sex_check.png",
        pedcheck="qc/peddy/{family}.ped_check.csv",
        pedcheck_png="qc/peddy/{family}.ped_check.png",
        peddy_ped="qc/peddy/{family}.peddy.ped",
        rel_diff="qc/peddy/{family}.ped_check.rel-difference.csv"
    log:
        "logs/qc/peddy/{family}.log"
    conda:
        "../envs/peddy.yaml"
    shell:
        '''
        mkdir -p qc/peddy
        peddy \
          --prefix ./qc/peddy/{wildcards.family} \
          --plot \
          {input.vcf} \
          {input.ped} \
          2>&1 | tee {log}
        '''

rule peddy_relatedness_mqc:
    input:
        pedcheck="qc/peddy/{family}.ped_check.csv"
    output:
        tsv="qc/multiqc_custom/{family}/peddy_relatedness_mqc.tsv"
    log:
        "logs/qc/peddy/{family}.relatedness_mqc.log"
    shell:
        '''
        mkdir -p qc/multiqc_custom/{wildcards.family}
        awk -F',' '
        BEGIN {{
            OFS="\t";
            print "Pair_ID","Sample_A","Sample_B","Peddy_Relatedness"
        }}
        NR>1 {{
            pid = $1 "_" $2
            rel = $3
            print pid, $1, $2, rel
        }}
        ' {input.pedcheck} > {output.tsv}
        '''

rule bcftools_stats:
    input:
        vcf=get_sequence_var_vcf
    output:
        stats="qc/bcftools/{family}/{sample}.stats"
    log: 
        "logs/qc/bcftools/{family}/{sample}.stats.log"
    conda:
        "../envs/common.yaml"
    shell:
        '''
        mkdir -p qc/bcftools/{wildcards.family} logs/qc/bcftools/{wildcards.family}
        bcftools stats \
            -s {wildcards.sample} \
            {input.vcf} \
        | awk -v sample="{wildcards.sample}" '
            BEGIN {{ OFS="\t" }}
            $1=="ID" && $2=="0" {{ $3=sample; print; next }}
            $1=="QUAL" {{ $3=int($3+0.5); print; next }}
            {{ print }}
            ' \
            > {output.stats}
        '''

rule verifybam:
    input:
        cram=get_cram
    output:
        selfsm="qc/verifybam/{family}/{sample}.selfSM"
    log:
        "logs/qc/verifybam/{family}/{sample}.verifybam.log"
    params:
        out_prefix="qc/verifybam/{family}/{sample}",
        sample="{sample}",
        ref=config["ref"]["genome"],
        svdp=config["qc"]["svdp"]
    wrapper:
        get_wrapper_path("verifybamid")

rule qualimap:
    input:
        bam="qc/qualimap/{family}/{sample}.bam"
    output:
        html="qc/qualimap/{family}/{sample}/qualimapReport.html",
        raw="qc/qualimap/{family}/{sample}/raw_data_qualimapReport/coverage_histogram.txt"
    log:
        "logs/qc/qualimap/{family}/{sample}.log"
    conda:
        "../envs/qualimap.yaml"
    threads: 8
    resources: 
        mem_mb=60000
    params:
        outdir="qc/qualimap/{family}/{sample}",
        java_mem="48G",
        nw=400,
        hm=3
    shell:
        '''
        mkdir -p {params.outdir}
        mkdir -p logs/qc/qualimap/{wildcards.family}
        mkdir -p {params.outdir}/tmp

        unset DISPLAY
        export _JAVA_OPTIONS="-Djava.io.tmpdir=$(pwd)/{params.outdir}/tmp"

        qualimap \
            --java-mem-size={params.java_mem} \
            bamqc \
            -bam {input.bam} \
            -outdir {params.outdir} \
            -nt {threads} \
            -c \
            -nw {params.nw} \
            -hm {params.hm} \
        &> {log}
        '''

rule dragenmetrics:
    input:
        samples_tsv=config["run"]["samples"],
        metrics_files=get_dragen_metrics_files
    output:
        duplication="qc/multiqc_custom/{family}/dragen_duplication_mqc.tsv",
        mapping="qc/multiqc_custom/{family}/dragen_mapping_mqc.tsv",
        qc_summary="qc/multiqc_custom/{family}/dragen_qc_summary_mqc.tsv"
    log:
        "logs/qc/dragenmetrics/{family}.log"
    conda:
        "../envs/annotate.yaml"
    params:
        outdir="qc/multiqc_custom/{family}"
    script:
        "../scripts/dragen_metrics_to_mqc.py"

rule multiqc:
    input:
        peddy_html=f"qc/peddy/{family}.html",
        peddy_relatedness="qc/multiqc_custom/{family}/peddy_relatedness_mqc.tsv",
        bcftools_stats=expand("qc/bcftools/{family}/{sample}.stats", family=family, sample=samples.index),
        selfsm=expand("qc/verifybam/{family}/{sample}.selfSM", family=family, sample=samples.index),
        qualimap_report=expand("qc/qualimap/{family}/{sample}/qualimapReport.html", family=family, sample=samples.index),
        dragen_duplication="qc/multiqc_custom/{family}/dragen_duplication_mqc.tsv",
        dragen_mapping="qc/multiqc_custom/{family}/dragen_mapping_mqc.tsv",
        dragen_qc_summary="qc/multiqc_custom/{family}/dragen_qc_summary_mqc.tsv",
    output:
        report="qc/multiqc/{family}.multiqc_report.html"
    log:
        "logs/qc/multiqc/{family}.multiqc.log"
    conda:
        "../envs/multiqc.yaml"
    shell:
        '''
        mkdir -p qc/multiqc
        multiqc qc \
        --force \
        --filename {wildcards.family}.multiqc_report.html \
        --config {workflow.basedir}/scripts/multiqc_config.yaml \
        -o qc/multiqc \
        &> {log}
        '''

rule publish_multiqc_report:
    input:
        report="qc/multiqc/{family}.multiqc_report.html"
    output:
        published="reports/{family}.multiqc_report.html"
    log:
        "logs/qc/multiqc/{family}.publish.log"
    shell:
        """
        mkdir -p reports
        cp {input.report} {output.published}
        """
