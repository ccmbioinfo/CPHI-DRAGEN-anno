rule snpeff:
    input: get_sv_vcf
    output:
        vcf = "sv/{family}.sv.snpeff.vcf",
    log:
        "logs/sv/{family}.snpeff.log"
    params:
        java_opts = config["params"]["snpeff"]["java_opts"],
        reference = config["ref"]["name"],
        data_dir = config["annotation"]["snpeff"]["data_dir"]
    wrapper:
        get_wrapper_path("snpeff")

rule annotsv:
    input: get_sv_vcf
    output:
        annotsv_annotated =  "sv/{family}.AnnotSV.tsv",
        annotsv_unannotated =  temp("sv/{family}.AnnotSV.unannotated.tsv")
    log: "logs/sv/{family}.annotsv.log"
    params:
        annotsv_path = config["tools"]["annotSV"]
    conda:
        "../envs/common.yaml"
    shell:
        """
        (export ANNOTSV={params.annotsv_path}; 
        {params.annotsv_path}/bin/AnnotSV \
            -SVinputFile {input} \
            -outputFile {output.annotsv_annotated} \
            -overlap 50 \
            -reciprocal 1 \
            -genomeBuild GRCh38) > {log} 2>&1
        """

rule sv_report:
    input: 
        snpeff = "sv/{family}.sv.snpeff.vcf",
        annotsv = "sv/{family}.AnnotSV.tsv"
    output: "sv/{family}.sv.csv"
    log: "logs/sv/{family}.sv.report.log"
    params:
        cphi_dragen = config["tools"]["cphi-dragen-anno"],
        HPO = config["run"]["hpo"] if config["run"]["hpo"] else "none",
        omim = config["annotation"]["omim_path"],
        exon = config["annotation"]["general"]["exon"],
        repeats = config["annotation"]["general"]["adotto_repeats"],
        gnomad_SV = config["annotation"]["sv_report"]["gnomad_SV"],
        dgv = config["annotation"]["sv_report"]["dgv"],
        ensembl = config["annotation"]["general"]["ensembl"],
        clingen_path = config["annotation"]["general"]["clingen_path"],
        samples = config["run"]["samples"]
    conda:
        "../envs/str_sv.yaml"
    shell:
        """
        if [[ {params.HPO} == "none" ]]
            then
                    (python3 {params.cphi_dragen}/workflow/scripts/annotate_SVs.py \
                        -annotsv {input.annotsv} \
                        -snpeff {input.snpeff} \
                        -variant_type SV \
                        -omim {params.omim} \
                        -exon {params.exon} \
                        -gnomad {params.gnomad_SV} \
                        -dgv {params.dgv} \
                        -ensembl {params.ensembl} \
                        -repeats {params.repeats} \
                        -clingen_HI {params.clingen_path}/ClinGen_haploinsufficiency_gene_GRCh38.bed \
                        -clingen_TS {params.clingen_path}/ClinGen_triplosensitivity_gene_GRCh38.bed \
                        -clingen_disease {params.clingen_path}/ClinGen_tableExport_202310.csv \
                        -clingen_regions {params.clingen_path}/ClinGen_region_curation_list_GRCh38.tsv \
                        -samples {params.samples}) > {log} 2>&1
            else
                (python3 {params.cphi_dragen}/workflow/scripts/annotate_SVs.py \
                    -annotsv {input.annotsv} \
                    -snpeff {input.snpeff} \
                    -variant_type SV \
                    -omim {params.omim} \
                    -hpo {params.HPO} \
                    -exon {params.exon} \
                    -gnomad {params.gnomad_SV} \
                    -dgv {params.dgv} \
                    -ensembl {params.ensembl} \
                    -repeats {params.repeats} \
                    -clingen_HI {params.clingen_path}/ClinGen_haploinsufficiency_gene_GRCh38.bed \
                    -clingen_TS {params.clingen_path}/ClinGen_triplosensitivity_gene_GRCh38.bed \
                    -clingen_disease {params.clingen_path}/ClinGen_tableExport_202310.csv \
                    -clingen_regions {params.clingen_path}/ClinGen_region_curation_list_GRCh38.tsv \
                    -samples {params.samples}) > {log} 2>&1
            fi
        """
