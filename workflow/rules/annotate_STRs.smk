rule repeat_VCF_to_df:
    input: get_repeat_dir
    output: "STRs/{family}.repeats.tsv"
    params:
        cphi_dragen_anno = config["tools"]["cphi-dragen-anno"]
    log: "logs/repeat_VCF_to_df/{family}.repeats.log"
    conda: "../envs/annotate.yaml"
    shell: "(python3 {params.cphi_dragen_anno}/workflow/scripts/repeat_VCF_to_df.py --repeat_vcf_dir {input} --output_file {output}) > {log} 2>&1"