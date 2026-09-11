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


rule repeat_VCF_to_df:
    input:
        samples_tsv=config["run"]["samples"],
        expansionhunter_vcfs=expand("STRs/expansionhunter/{sample}.vcf", sample=samples.index),
        variant_catalog=config["annotation"]["str_variant_catalog"],
        disease_thresholds=config["annotation"]["str_disease_thresholds"]
    output: temp("STRs/{family}.repeats.tsv")
    params:
        cphi_dragen_anno=config["tools"]["cphi-dragen-anno"],
        expansionhunter_dir="STRs/expansionhunter"
    log: "logs/STRs/{family}.repeats.log"
    conda: "../envs/annotate.yaml"
    shell:
        """
        python3 {params.cphi_dragen_anno}/workflow/scripts/repeat_VCF_to_df.py \
            --samples_tsv {input.samples_tsv} \
            --family {wildcards.family} \
            --expansionhunter_dir {params.expansionhunter_dir} \
            --variant_catalog {input.variant_catalog} \
            --disease_thresholds {input.disease_thresholds} \
            --output_file {output} \
            > {log} 2>&1
        """

rule annotate_path_str_loci:
    input: 
        repeat_tsv="STRs/{family}.repeats.tsv",
        samples_tsv=config["run"]["samples"]
    output: "reports/{family}.known.path.str.loci.hg38.csv"
    params:
        cphi_dragen_anno = config["tools"]["cphi-dragen-anno"],
        disease_thresholds = config["annotation"]["str_disease_thresholds"]
    log: "logs/STRs/{family}.STR.report.log"
    conda: "../envs/annotate.yaml"
    shell: "(python3 {params.cphi_dragen_anno}/workflow/scripts/annotate_path_str_loci.py --repeat_tsv {input.repeat_tsv} --disease_thresholds {params.disease_thresholds} --samples_tsv {input.samples_tsv} --output_file {output}) > {log} 2>&1"
