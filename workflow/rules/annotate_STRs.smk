rule repeat_VCF_to_df:
    input: config["run"]["samples"]
    output: "STRs/{family}.repeats.tsv"
    params:
        cphi_dragen_anno = config["tools"]["cphi-dragen-anno"]
    log: "logs/STRs/{family}.repeats.log"
    conda: "../envs/annotate.yaml"
    shell: "(python3 {params.cphi_dragen_anno}/workflow/scripts/repeat_VCF_to_df.py --samples_tsv {input} --output_file {output}) > {log} 2>&1"

rule annotate_path_str_loci:
    input: "STRs/{family}.repeats.tsv"
    output: "reports/{family}.known.path.str.loci.hg38.csv"
    params:
        cphi_dragen_anno = config["tools"]["cphi-dragen-anno"],
        disease_thresholds = config["annotation"]["str_disease_thresholds"]
    log: "logs/STRs/{family}.STR.report.log"
    conda: "../envs/annotate.yaml"
    shell: "(python3 {params.cphi_dragen_anno}/workflow/scripts/annotate_path_str_loci.py --repeat_tsv {input} --disease_thresholds {params.disease_thresholds} --output_file {output}) > {log} 2>&1" 