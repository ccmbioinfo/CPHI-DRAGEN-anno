rule merge_VCFs:
    input: get_repeat_dir
    output: "STRs/{family}.repeats.vcf.gz"
    log: "logs/STRs/{family}.repeats.log"
    conda: "../envs/common.yaml"
    shell: "(bcftools merge -O z -o {output} {input}/*.repeats.vcf.gz) > {log} 2>&1"

rule repeat_VCF_to_df:
    input: 
        vcf = "STRs/{family}.repeats.vcf.gz",
        tbi = "STRs/{family}.repeats.vcf.gz.tbi"
    output: "STRs/{family}.repeats.tsv"
    params:
        cphi_dragen_anno = config["tools"]["cphi-dragen-anno"]
    log: "logs/repeat_VCF_to_df/{family}.repeats.log"
    conda: "../envs/annotate.yaml"
    shell: "(python3 {params.cphi_dragen_anno}/workflow/scripts/repeat_VCF_to_df.py --repeat_vcf {input.vcf} --output_file {output}) > {log} 2>&1"

rule tabix:
    input: 
        "{prefix}.vcf.gz"
    output: 
        "{prefix}.vcf.gz.tbi"
    log: 
        "logs/{prefix}.log"
    conda:
           "../envs/common.yaml"
    shell: "(tabix {input}) > {log} 2>&1"