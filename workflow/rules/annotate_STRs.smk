rule expansionhunter:
    input:
        cram=get_cram,
        reference=config["ref"]["genome"],
        catalog=config["annotation"]["str_variant_catalog"],
        sex_check=f"qc/peddy/{family}.sex_check.csv"
    output:
        vcf=temp("STRs/expansionhunter/{sample}.vcf"),
        json=temp("STRs/expansionhunter/{sample}.json"),
        realigned_bam=temp("STRs/expansionhunter/{sample}_realigned.bam")
    params:
        output_prefix="STRs/expansionhunter/{sample}",
        sif=config["tools"]["expansionhunter_sif"]
    log:
        "logs/STRs/expansionhunter/{sample}.log"
    conda:
        "../envs/expansionhunter.yaml"
    shell:
        """
        mkdir -p STRs/expansionhunter logs/STRs/expansionhunter
        sex=$(python3 -c 'import csv, sys; print(next(row["ped_sex"] for row in csv.DictReader(open(sys.argv[1])) if row["sample_id"] == sys.argv[2]))' {input.sex_check} {wildcards.sample})
        sex_arg=""
        if [[ "$sex" == "male" || "$sex" == "female" ]]; then
            sex_arg="--sex $sex"
        fi
        work_dir=$(pwd)
        cram_dir=$(dirname {input.cram})
        reference_dir=$(dirname {input.reference})
        catalog_dir=$(dirname {input.catalog})
        echo "ExpansionHunter ped_sex for {wildcards.sample}: $sex" > {log}
        apptainer exec \
            --bind "$work_dir:$work_dir" \
            --bind "$cram_dir:$cram_dir:ro" \
            --bind "$reference_dir:$reference_dir:ro" \
            --bind "$catalog_dir:$catalog_dir:ro" \
            --pwd "$work_dir" \
            {params.sif} ExpansionHunter \
            --reads {input.cram} \
            --reference {input.reference} \
            --variant-catalog {input.catalog} \
            --output-prefix {params.output_prefix} \
            $sex_arg \
            >> {log} 2>&1
        """


rule sort_expansionhunter_bam:
    input:
        bam="STRs/expansionhunter/{sample}_realigned.bam"
    output:
        bam=temp("STRs/expansionhunter/{sample}_realigned.sorted.bam"),
        bai=temp("STRs/expansionhunter/{sample}_realigned.sorted.bam.bai")
    log:
        "logs/STRs/reviewer/{sample}.sort.log"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        mkdir -p logs/STRs/reviewer
        samtools sort -o {output.bam} {input.bam} > {log} 2>&1 && samtools index {output.bam} {output.bai} >> {log} 2>&1
        """


rule reviewer:
    input:
        bam="STRs/expansionhunter/{sample}_realigned.sorted.bam",
        bai="STRs/expansionhunter/{sample}_realigned.sorted.bam.bai",
        vcf="STRs/expansionhunter/{sample}.vcf",
        reference=config["ref"]["genome"],
        reference_fai=config["ref"]["genome"] + ".fai",
        catalog=config["annotation"]["str_variant_catalog"]
    output:
        reviewer_dir=temp(directory("STRs/reviewer_work/{sample}"))
    params:
        log_dir="logs/STRs/reviewer/{sample}"
    conda:
        "../envs/str_tools.yaml"
    shell:
        """
        mkdir -p {output.reviewer_dir} {params.log_dir}
        for locus in $(python3 -c 'import json, sys; print(*(row["LocusId"] for row in json.load(open(sys.argv[1]))), sep="\\n")' {input.catalog}); do
            mkdir -p {output.reviewer_dir}/$locus {params.log_dir}/$locus
            REViewer \
                --reads {input.bam} \
                --vcf {input.vcf} \
                --reference {input.reference} \
                --catalog {input.catalog} \
                --locus $locus \
                --output-prefix {output.reviewer_dir}/$locus/{wildcards.sample} \
                > {params.log_dir}/$locus/{wildcards.sample}.reviewer.log 2>&1 || true
        done
        """


rule repeat_VCF_to_df:
    input:
        samples_tsv=config["run"]["samples"],
        expansionhunter_vcfs=expand("STRs/expansionhunter/{sample}.vcf", sample=samples.index),
        disease_thresholds=config["annotation"]["str_disease_thresholds"]
    output:
        temp("STRs/{family}.repeats.tsv")
    params:
        cphi_dragen_anno=config["tools"]["cphi-dragen-anno"],
        expansionhunter_dir="STRs/expansionhunter"
    log:
        "logs/STRs/{family}.repeats.log"
    conda:
        "../envs/annotate.yaml"
    shell:
        """
        python3 {params.cphi_dragen_anno}/workflow/scripts/repeat_VCF_to_df.py \
            --samples_tsv {input.samples_tsv} \
            --expansionhunter_dir {params.expansionhunter_dir} \
            --disease_thresholds {input.disease_thresholds} \
            --output_file {output} \
            > {log} 2>&1
        """


rule annotate_path_str_loci:
    input:
        repeat_tsv="STRs/{family}.repeats.tsv",
        samples_tsv=config["run"]["samples"],
        disease_thresholds=config["annotation"]["str_disease_thresholds"]
    output:
        "reports/{family}.known.path.str.loci.hg38.csv"
    params:
        cphi_dragen_anno=config["tools"]["cphi-dragen-anno"]
    log:
        "logs/STRs/{family}.STR.report.log"
    conda:
        "../envs/annotate.yaml"
    shell:
        """
        python3 {params.cphi_dragen_anno}/workflow/scripts/annotate_path_str_loci.py \
            --repeat_tsv {input.repeat_tsv} \
            --disease_thresholds {input.disease_thresholds} \
            --samples_tsv {input.samples_tsv} \
            --output_file {output} \
            > {log} 2>&1
        """


rule build_str_reviewer_site:
    input:
        report="reports/{family}.known.path.str.loci.hg38.csv",
        samples_tsv=config["run"]["samples"],
        catalog=config["annotation"]["str_variant_catalog"],
        reviewer_dirs=expand("STRs/reviewer_work/{sample}", sample=samples.index)
    output:
        metadata="STRs/reviewer/{family}/flipbook_metadata.tsv",
        site=directory("STRs/reviewer/{family}/site"),
        launcher="reports/{family}.STR.review.html"
    log:
        "logs/STRs/{family}.reviewer_site.log"
    conda:
        "../envs/str_tools.yaml"
    script:
        "../scripts/str/build_str_reviewer_site.py"
