
crg2_pacbio = config["tools"]["crg2_pacbio"]

rule extract_mt_variants:
    input:
        vcf=get_small_variant_vcf
    output:
        temp("mitochondrial_variants/{family}.dragen.mt.vcf.gz")
    params:
        outdir="mitochondrial_variants/",
        mt_contig="chrM"
    log:
        "logs/mity/extract_mt_from_dragen_vcf/{family}.extract_mt_from_dragen_vcf.log"
    wrapper:    
        get_wrapper_path("mito/extract_from_dragen_vcf")

rule bcftools_normalise:
    input:
        "mitochondrial_variants/{family}.dragen.mt.vcf.gz"
    output:
        temp("mitochondrial_variants/{family}.mt.normalise.decompose.vcf.gz")
    params:
        outdir="mitochondrial_variants/",
        tool=config["tools"]["mity"],
        cphi_dragen_anno=config["tools"]["cphi-dragen-anno"]
    log:
        "logs/mity/bcftools_normalise/{family}.bcftools_normalise.log"
    wrapper:    
        get_wrapper_path("mito/bcftools_normalise")

rule mity_report:
    input:
        "mitochondrial_variants/{family}.mt.normalise.decompose.vcf.gz"
    output:
        "mitochondrial_variants/{family}.mity.annotated.vcf.gz",
        "mitochondrial_variants/{family}.mity.report.xlsx"
    params:
        outdir="mitochondrial_variants/",
        tool=config["tools"]["mity"],
        report_config=config["annotation"]["mity"]["report_config"],
        vcfanno_config=config["annotation"]["mity"]["vcfanno_config"],
        base_path=config["annotation"]["mity"]["base_path"]
    log:
        "logs/mity/mity_report/{family}.mity_report.log"
    wrapper:
        get_wrapper_path("mito/report")

rule generate_mt_report:
    input:
        vcf="mitochondrial_variants/{family}.mity.annotated.vcf.gz",
        report="mitochondrial_variants/{family}.mity.report.xlsx"
    output:
        "reports/{family}.mito.hg38.csv"        
    log:
        "logs/report/mitochondrial/{family}.mitochondrial.report.log"
    conda:
        "../envs/mt_report.yaml"
    params:
        crg2_pacbio=crg2_pacbio,
        vaf_field="AF",
        report_suffix=".hg38"
    script:
        crg2_pacbio + "/scripts/mt_report.py"
