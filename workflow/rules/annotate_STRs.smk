rule repeat_VCF_to_df:
    input: config["run"]["samples"]
    output: "STRs/{family}.repeats.tsv"
    params:
        cphi_dragen_anno = config["tools"]["cphi-dragen-anno"]
    log: "logs/STRs/{family}.repeats.log"
    conda: "../envs/annotate.yaml"
    shell: "(python3 {params.cphi_dragen_anno}/workflow/scripts/repeat_VCF_to_df.py --samples_tsv {input} --family {wildcards.family} --output_file {output}) > {log} 2>&1"

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