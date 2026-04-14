rule get_sequence_variants_for_CH:
    input:
        gemini_db="annotated/coding/{family}-gemini.db"
    output:
        variants="sequence_variants/{family}.{severity}.impact.variants.tsv",
    params:
        severity="{severity}",
        crg2_pacbio = config["tools"]["crg2_pacbio"],
        seq_type="short"
    log:
        "logs/compound_hets/{family}.get.sequence.variants.for.CH.{severity}.log",
    conda:
        "../envs/gemini.yaml"
    wildcard_constraints:
        severity="HIGH-MED|LOW"
    shell:
        "{params.crg2_pacbio}/scripts/compound_hets/get_sequence_var_for_CH.sh {input.gemini_db} {params.severity} {params.seq_type} > {output.variants}"

rule get_VCF_sample_order:
    input:
        vcf="annotated/coding/vcfanno/{family}.coding.vep.vcfanno.vcf.gz",
    output:
        sample_order="sequence_variants/{family}.sample.order.txt",
    log:
        "logs/compound_hets/{family}.get.VCF.sample.order.log",
    conda:
        "../envs/common.yaml"
    shell:
        "bcftools query -l {input.vcf} > {output.sample_order}"

if len(children) > 0:
    rule identify_compound_hets_with_denovo:
        input:
            high_med_variants="sequence_variants/{family}.HIGH-MED.impact.variants.tsv",
            low_variants="sequence_variants/{family}.LOW.impact.variants.tsv",
            small_variant_report_dir="sequence_variants/coding/{family}",
            panel_variant_report_dir="sequence_variants/panel/{family}",
            panel_flank_variant_report_dir="sequence_variants/panel-flank/{family}",
            wgs_high_impact_variant_report_dir="sequence_variants/wgs-high-impact/{family}",
            wgs_denovo_variant_report_dir="sequence_variants/denovo/{family}",
            SV_report="sv/{family}.sv.csv",
            CNV_report="cnv/{family}.cnv.csv",
            ensembl=config["annotation"]["general"]["ensembl"],
            ensembl_to_NCBI_df=config["annotation"]["general"]["ensembl_to_NCBI_df"],
            HPO=config["run"]["hpo"],
            pedigree=config["run"]["ped"],
            sample_order="sequence_variants/{family}.sample.order.txt",
        output:
            sequence_variant_report_CH="reports/{family}.wgs.coding.CH.csv",
            panel_variant_report_CH="reports/{family}.panel.CH.csv",
            panel_flank_variant_report_CH="reports/{family}.panel-flank.CH.csv",
            wgs_high_impact_variant_report_CH="reports/{family}.wgs.high.impact.CH.csv",
            wgs_denovo_variant_report_CH="reports/{family}.wgs.denovo.CH.csv",
            SV_report_CH="reports/{family}.sv.CH.csv",
            CNV_report_CH="reports/{family}.cnv.CH.csv",
            compound_het_status="reports/{family}.compound.het.status.CH.csv",
        params:
            crg2_pacbio = config["tools"]["crg2_pacbio"],
            seq_type="short"
        conda:
            "../envs/str_sv.yaml"
        log:
            "logs/compound_hets/{family}.identify.compound.hets.log",
        shell:
            """
            (python3 {params.crg2_pacbio}/scripts/annotate_compound_hets.py --seq_type {params.seq_type} --high_med {input.high_med_variants} \
            --low {input.low_variants} \
            --sv {input.SV_report}  \
            --cnv {input.CNV_report}  \
            --ensembl {input.ensembl}  \
            --ensembl_to_NCBI_df {input.ensembl_to_NCBI_df}  \
            --pedigree {input.pedigree}  \
            --hpo {input.HPO}  \
            --sequence_variant_report_dir {input.small_variant_report_dir}  \
            --panel_variant_report_dir {input.panel_variant_report_dir}  \
            --panel_flank_variant_report_dir {input.panel_flank_variant_report_dir}  \
            --wgs_high_impact_variant_report_dir {input.wgs_high_impact_variant_report_dir}  \
            --wgs_denovo_variant_report_dir {input.wgs_denovo_variant_report_dir}  \
            --sample_order {input.sample_order}  \
            --family {wildcards.family}) > {log} 2>&1
            """
else:
        rule identify_compound_hets:
            input:
                high_med_variants="sequence_variants/{family}.HIGH-MED.impact.variants.tsv",
                low_variants="sequence_variants/{family}.LOW.impact.variants.tsv",
                small_variant_report_dir="sequence_variants/coding/{family}",
                panel_variant_report_dir="sequence_variants/panel/{family}",
                panel_flank_variant_report_dir="sequence_variants/panel-flank/{family}",
                wgs_high_impact_variant_report_dir="sequence_variants/wgs-high-impact/{family}",
                SV_report="reports/{family}.sv.csv",
                CNV_report="reports/{family}.cnv.csv",
                ensembl=config["annotation"]["general"]["ensembl"],
                ensembl_to_NCBI_df=config["annotation"]["general"]["ensembl_to_NCBI_df"],
                HPO=config["run"]["hpo"],
                pedigree=config["run"]["ped"],
                sample_order="sequence_variants/{family}.sample.order.txt",
            output:
                sequence_variant_report_CH="reports/{family}.wgs.coding.CH.csv",
                panel_variant_report_CH="reports/{family}.panel.CH.csv",
                panel_flank_variant_report_CH="reports/{family}.panel-flank.CH.csv",
                wgs_high_impact_variant_report_CH="reports/{family}.wgs.high.impact.CH.csv",
                SV_report_CH="reports/{family}.sv.CH.csv",
                CNV_report_CH="reports/{family}.cnv.CH.csv",
                compound_het_status="reports/{family}.compound.het.status.CH.csv",
            params:
                crg2_pacbio = config["tools"]["crg2_pacbio"],
                seq_type="short"
            conda:
                "../envs/str_sv.yaml"
            log:
                "logs/compound_hets/{family}.identify.compound.hets.log",
            shell:
                """
                (python3 {params.crg2_pacbio}/scripts/annotate_compound_hets.py --seq_type {params.seq_type} --high_med {input.high_med_variants} \
                --low {input.low_variants} \
                --sv {input.SV_report}  \
                --cnv {input.CNV_report}  \
                --ensembl {input.ensembl}  \
                --ensembl_to_NCBI_df {input.ensembl_to_NCBI_df}  \
                --pedigree {input.pedigree}  \
                --hpo {input.HPO}  \
                --sequence_variant_report_dir {input.small_variant_report_dir}  \
                --panel_variant_report_dir {input.panel_variant_report_dir}  \
                --panel_flank_variant_report_dir {input.panel_flank_variant_report_dir}  \
                --wgs_high_impact_variant_report_dir {input.wgs_high_impact_variant_report_dir}  \
                --sample_order {input.sample_order}  \
                --family {wildcards.family}) > {log} 2>&1
                """
 