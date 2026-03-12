rule cram_to_bam:
    input:
        cram=get_cram
    output:
        bam=temp("qc/picard/{family}/{sample}.bam")
    log:
        "logs/qc/picard/{family}/{sample}.cram_to_bam.log"
    conda:
        "../envs/samtools.yaml"
    threads: 8
    params:
        ref=config["ref"]["genome"]
    shell:
        """
        mkdir -p qc/picard/{wildcards.family} logs/qc/picard/{wildcards.family}
        samtools view -@ {threads} -T {params.ref} -b {input.cram} > {output.bam}
        """

rule picard_markduplicates:
    input:
        bam="qc/picard/{family}/{sample}.bam"
    output:
        bam=temp("qc/picard/{family}/{sample}.markdup.bam"),
        metrics="qc/picard/{family}/{sample}.duplication_metrics.txt"
    log:
        "logs/qc/picard/{family}/{sample}.markdup.log"
    conda:
        "../envs/picard.yaml"
    threads: 4
    resources:
        mem_mb=60000
    params:
        java_mem="48G",
        tmpdir="qc/picard/{family}/{sample}.tmp"
    shell:
        """
        mkdir -p {params.tmpdir}
        picard -Xmx{params.java_mem} MarkDuplicates \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.metrics} \
            TMP_DIR={params.tmpdir} \
            ASSUME_SORT_ORDER=coordinate \
            VALIDATION_STRINGENCY=SILENT \
            CREATE_INDEX=false \
            &> {log}
        """

rule peddy:
    input:
        vcf=get_smallvariants_vcf,
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
        vcf=get_smallvariants_vcf
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

rule samtools_stats:
    input:
        cram=get_cram
    output:
        stats="qc/samtools/{family}/{sample}.stats"
    log:
        "logs/qc/samtools/{family}/{sample}.log"
    conda:
        "../envs/samtools.yaml"
    params:
        ref=config["ref"]["genome"],
        refcache=config["qc"]["refcache"]
    shell:
        '''
        mkdir -p qc/samtools/{wildcards.family} logs/qc/samtools/{wildcards.family}
        export REF_CACHE="{params.refcache}/%2s/%2s/%s"
        samtools stats -r {params.ref} {input.cram} > {output.stats} 
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
        bam="qc/picard/{family}/{sample}.bam"
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

rule multiqc:
    input:
        peddy_html=f"qc/peddy/{family}.html",
        peddy_relatedness="qc/multiqc_custom/{family}/peddy_relatedness_mqc.tsv",
        bcftools_stats=expand("qc/bcftools/{family}/{sample}.stats", family=family, sample=samples.index),
        selfsm=expand("qc/verifybam/{family}/{sample}.selfSM", family=family, sample=samples.index),
        samtools_stats=expand("qc/samtools/{family}/{sample}.stats", family=family, sample=samples.index),
        qualimap_report=expand("qc/qualimap/{family}/{sample}/qualimapReport.html", family=family, sample=samples.index),
        picard_report=expand("qc/picard/{family}/{sample}.duplication_metrics.txt", family=family, sample=samples.index)
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
