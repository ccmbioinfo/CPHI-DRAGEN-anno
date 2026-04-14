rule cnv_snpeff:
    input: get_cnv_vcf
    output:
        vcf = "cnv/{family}.cnv.snpeff.vcf",
    log:
        "logs/cnv/{family}.snpeff.log"
    params:
        java_opts = config["params"]["snpeff"]["java_opts"],
        reference = config["ref"]["name"],
        data_dir = config["annotation"]["snpeff"]["data_dir"],
        config_file = config["annotation"]["snpeff"]["config"]
    wrapper:
        get_wrapper_path("snpeff")

rule cnv_annotsv:
    input: get_cnv_vcf
    output:
        annotsv_annotated =  "cnv/{family}.AnnotSV.tsv",
        annotsv_unannotated =  temp("cnv/{family}.AnnotSV.unannotated.tsv")
    log: "logs/cnv/{family}.annotsv.log"
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

rule cnv_report:
    input: 
        snpeff = "cnv/{family}.cnv.snpeff.vcf",
        annotsv = "cnv/{family}.AnnotSV.tsv"
    output: "cnv/{family}.cnv.csv"
    log: "logs/cnv/{family}.cnv.report.log"
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
        samples = config["run"]["samples"],
        thousandg = config["annotation"]["sv_report"]["1000G_SV"]
    conda:
        "../envs/str_sv.yaml"
    shell:
        """
        if [[ {params.HPO} == "none" ]]
            then
                    (python3 {params.cphi_dragen}/workflow/scripts/annotate_SVs.py \
                        -annotsv {input.annotsv} \
                        -snpeff {input.snpeff} \
                        -variant_type CNV \
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
                        -thousandg {params.thousandg} \
                        -samples {params.samples}) > {log} 2>&1

            else
                (python3 {params.cphi_dragen}/workflow/scripts/annotate_SVs.py \
                    -annotsv {input.annotsv} \
                    -snpeff {input.snpeff} \
                    -variant_type CNV \
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
                    -thousandg {params.thousandg} \
                    -samples {params.samples}) > {log} 2>&1
            fi
        """