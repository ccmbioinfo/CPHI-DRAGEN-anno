def output_status(output_path):
    if str(config["run"].get("acmg_sf", "")).lower() == "true":
        return temp(output_path)
    return(output_path)

hpo_available = config["run"].get("hpo", "")

hpo_panel_inputs = {
    "panel_variant_report_dir": "small_variants/panel/{family}",
    "panel_flank_variant_report_dir": "small_variants/panel-flank/{family}",
    "HPO": config["run"]["hpo"],
} if hpo_available else {}

hpo_panel_outputs = {
    "panel_variant_report_CH": "reports/{family}.panel.CH.hg38.csv",
    "panel_flank_variant_report_CH": "reports/{family}.panel-flank.CH.hg38.csv",
} if hpo_available else {}

def get_hpo_panel_args(wildcards, input):
    if hpo_available:
        return (
            f"--hpo {input.HPO} "
            f"--panel_variant_report_dir {input.panel_variant_report_dir} "
            f"--panel_flank_variant_report_dir {input.panel_flank_variant_report_dir}"
        )
    return ""

rule get_small_variants_for_CH:
    input:
        gemini_db="annotated/coding/{family}-gemini.db"
    output:
        variants=temp("small_variants/{family}.{severity}.impact.variants.tsv"),
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
        sample_order=temp("small_variants/{family}.sample.order.txt"),
    log:
        "logs/compound_hets/{family}.get.VCF.sample.order.log",
    conda:
        "../envs/common.yaml"
    shell:
        "bcftools query -l {input.vcf} > {output.sample_order}"

if len(children) > 0:
    rule identify_compound_hets_with_denovo:
        input:
            high_med_variants="small_variants/{family}.HIGH-MED.impact.variants.tsv",
            low_variants="small_variants/{family}.LOW.impact.variants.tsv",
            small_variant_report_dir="small_variants/coding/{family}",
            wgs_high_impact_variant_report_dir="small_variants/wgs-high-impact/{family}",
            wgs_denovo_variant_report_dir="small_variants/denovo/{family}",
            SV_report="sv/{family}.sv.csv",
            CNV_report="cnv/{family}.cnv.csv",
            ensembl=config["annotation"]["general"]["ensembl"],
            ensembl_to_NCBI_df=config["annotation"]["general"]["ensembl_to_NCBI_df"],
            pedigree=config["run"]["ped"],
            sample_order="small_variants/{family}.sample.order.txt",
            **hpo_panel_inputs,
        output:
            small_variant_report_CH=output_status("reports/{family}.wgs.coding.CH.hg38.csv"),
            wgs_high_impact_variant_report_CH=output_status("reports/{family}.wgs.high.impact.CH.hg38.csv"),
            wgs_denovo_variant_report_CH=output_status("reports/{family}.wgs.denovo.CH.hg38.csv"),
            SV_report_CH=output_status("reports/{family}.sv.CH.hg38.csv"),
            CNV_report_CH=output_status("reports/{family}.cnv.CH.hg38.csv"),
            compound_het_status="reports/{family}.compound.het.status.CH.hg38.csv",
            **hpo_panel_outputs,
        params:
            crg2_pacbio = config["tools"]["crg2_pacbio"],
            seq_type="short",
            hpo_panel_args=get_hpo_panel_args,
            acmg_sf_flag = str(config["run"].get("acmg_sf", "false")).lower(),
            mavedb_tsv = config["annotation"]["general"]["mavedb_tsv"]
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
            {params.hpo_panel_args}  \
            --sequence_variant_report_dir {input.small_variant_report_dir}  \
            --wgs_high_impact_variant_report_dir {input.wgs_high_impact_variant_report_dir}  \
            --wgs_denovo_variant_report_dir {input.wgs_denovo_variant_report_dir}  \
            --sample_order {input.sample_order}  \
            --family {wildcards.family}  \
            --acmg_sf {params.acmg_sf_flag} && \
            python3 {params.crg2_pacbio}/scripts/add_mavedb_columns.py \
            --family {wildcards.family} \
            --reports-dir reports \
            --mavedb-tsv {params.mavedb_tsv}) > {log} 2>&1
            """
else:
        rule identify_compound_hets:
            input:
                high_med_variants="small_variants/{family}.HIGH-MED.impact.variants.tsv",
                low_variants="small_variants/{family}.LOW.impact.variants.tsv",
                small_variant_report_dir="small_variants/coding/{family}",
                wgs_high_impact_variant_report_dir="small_variants/wgs-high-impact/{family}",
                SV_report="sv/{family}.sv.csv",
                CNV_report="cnv/{family}.cnv.csv",
                ensembl=config["annotation"]["general"]["ensembl"],
                ensembl_to_NCBI_df=config["annotation"]["general"]["ensembl_to_NCBI_df"],
                pedigree=config["run"]["ped"],
                sample_order="small_variants/{family}.sample.order.txt",
                **hpo_panel_inputs,
            output:
                small_variant_report_CH=output_status("reports/{family}.wgs.coding.CH.hg38.csv"),
                wgs_high_impact_variant_report_CH=output_status("reports/{family}.wgs.high.impact.CH.hg38.csv"),
                SV_report_CH=output_status("reports/{family}.sv.CH.hg38.csv"),
                CNV_report_CH=output_status("reports/{family}.cnv.CH.hg38.csv"),
                compound_het_status="reports/{family}.compound.het.status.CH.hg38.csv",
                **hpo_panel_outputs,
            params:
                crg2_pacbio = config["tools"]["crg2_pacbio"],
                seq_type="short",
                hpo_panel_args=get_hpo_panel_args,
                acmg_sf_flag = str(config["run"].get("acmg_sf", "false")).lower(),
                mavedb_tsv = config["annotation"]["general"]["mavedb_tsv"]
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
                {params.hpo_panel_args}  \
                --sequence_variant_report_dir {input.small_variant_report_dir}  \
                --wgs_high_impact_variant_report_dir {input.wgs_high_impact_variant_report_dir}  \
                --sample_order {input.sample_order}  \
                --family {wildcards.family}  \
                --acmg_sf {params.acmg_sf_flag} && \
                python3 {params.crg2_pacbio}/scripts/add_mavedb_columns.py \
                --family {wildcards.family} \
                --reports-dir reports \
                --mavedb-tsv {params.mavedb_tsv}) > {log} 2>&1
                """

slivar_hpo_panel_inputs = {
    "panel_variant_report": "small_variants_slivar/panel/{family}/{family}.panel.slivar.hg38.csv",
    "panel_flank_variant_report": "small_variants_slivar/panel-flank/{family}/{family}.panel-flank.slivar.hg38.csv",
    "HPO": config["run"]["hpo"],
} if hpo_available else {}

slivar_hpo_panel_outputs = {
    "panel_variant_report_CH": "reports_slivar/{family}.panel.CH.hg38.csv",
    "panel_flank_variant_report_CH": "reports_slivar/{family}.panel-flank.CH.hg38.csv",
} if hpo_available else {}

def get_slivar_hpo_panel_args(wildcards, input):
    if hpo_available:
        return (
            '--hpo "$hpo" '
            "--panel_variant_report_dir input/panel "
            "--panel_flank_variant_report_dir input/panel-flank"
        )
    return ""

def stage_slivar_hpo_panel_reports(wildcards, input):
    if hpo_available:
        return (
            f'hpo=$(realpath {input.HPO})\n'
            'mkdir -p "$stage/input/panel" "$stage/input/panel-flank"\n'
            f'ln -sf "$(realpath {input.panel_variant_report})" "$stage/input/panel/{wildcards.family}.panel.wgs.slivar.hg38.csv"\n'
            f'ln -sf "$(realpath {input.panel_flank_variant_report})" "$stage/input/panel-flank/{wildcards.family}.panel-flank.wgs.slivar.hg38.csv"'
        )
    return ""

rule slivar_select_compound_het_candidates:
    input:
        vcf="annotated/coding/vcfanno/{family}.coding.vep.vcfanno.vcf.gz"
    output:
        candidates=temp("small_variants_slivar/compound-hets/{family}/{family}.compound-het.candidates.vcf.gz"),
    log:
        "logs/compound_hets/{family}.slivar.select.compound.het.candidates.log",
    conda:
        "../wrappers/slivar/environment.yaml"
    params:
        js=f"{workflow.basedir}/scripts/slivar/slivar_functions.js",
        consequence_order_file=f"{workflow.basedir}/scripts/slivar/default-order.txt",
        mode="compound-hets",
    wrapper:
        get_wrapper_path("slivar")

rule get_slivar_small_variants_for_CH:
    input:
        vcf="small_variants_slivar/compound-hets/{family}/{family}.compound-het.candidates.vcf.gz",
    output:
        high_med=temp("small_variants_slivar/{family}.HIGH-MED.impact.variants.tsv"),
        low=temp("small_variants_slivar/{family}.LOW.impact.variants.tsv"),
    log:
        "logs/compound_hets/{family}.slivar.get.sequence.variants.for.CH.log",
    conda:
        "../envs/slivar.yaml"
    shell:
        """
        (set -e
        mkdir -p $(dirname {output.high_med})
        python3 {workflow.basedir}/scripts/slivar/build_ch_tsv.py \
        --vcf {input.vcf} \
        --impact-order-file {workflow.basedir}/scripts/slivar/default-order.txt \
        --high-med-out {output.high_med} \
        --low-out {output.low}) > {log} 2>&1
        """

if len(children) > 0:
    rule identify_compound_hets_slivar_with_denovo:
        input:
            high_med_variants="small_variants_slivar/{family}.HIGH-MED.impact.variants.tsv",
            low_variants="small_variants_slivar/{family}.LOW.impact.variants.tsv",
            small_variant_report="small_variants_slivar/coding/{family}/{family}.coding.slivar.hg38.csv",
            wgs_high_impact_variant_report="small_variants_slivar/wgs-high-impact/{family}/{family}.wgs-high-impact.slivar.hg38.csv",
            wgs_denovo_variant_report="small_variants_slivar/denovo/{family}/{family}.denovo.slivar.hg38.csv",
            SV_report="sv/{family}.sv.csv",
            CNV_report="cnv/{family}.cnv.csv",
            ensembl=config["annotation"]["general"]["ensembl"],
            ensembl_to_NCBI_df=config["annotation"]["general"]["ensembl_to_NCBI_df"],
            pedigree=config["run"]["ped"],
            sample_order="small_variants/{family}.sample.order.txt",
            **slivar_hpo_panel_inputs,
        output:
            small_variant_report_CH=output_status("reports_slivar/{family}.wgs.coding.CH.hg38.csv"),
            wgs_high_impact_variant_report_CH=output_status("reports_slivar/{family}.wgs.high.impact.CH.hg38.csv"),
            wgs_denovo_variant_report_CH=output_status("reports_slivar/{family}.wgs.denovo.CH.hg38.csv"),
            SV_report_CH=output_status("reports_slivar/{family}.sv.CH.hg38.csv"),
            CNV_report_CH=output_status("reports_slivar/{family}.cnv.CH.hg38.csv"),
            compound_het_status="reports_slivar/{family}.compound.het.status.CH.hg38.csv",
            **slivar_hpo_panel_outputs,
        params:
            crg2_pacbio=config["tools"]["crg2_pacbio"],
            seq_type="short",
            hpo_panel_args=get_slivar_hpo_panel_args,
            stage_hpo_panel_reports=stage_slivar_hpo_panel_reports,
            acmg_sf_flag=str(config["run"].get("acmg_sf", "false")).lower(),
            mavedb_tsv=config["annotation"]["general"]["mavedb_tsv"]
        conda:
            "../envs/str_sv.yaml"
        log:
            "logs/compound_hets/{family}.slivar.identify.compound.hets.log",
        shell:
            """
            (set -e
            root=$(pwd)
            stage=$(mktemp -d "$root/.slivar_ch_work.{wildcards.family}.XXXXXX")
            trap 'rm -rf "$stage"' EXIT
            high_med=$(realpath {input.high_med_variants})
            low=$(realpath {input.low_variants})
            sample_order=$(realpath {input.sample_order})
            ensembl=$(realpath {input.ensembl})
            ensembl_to_ncbi=$(realpath {input.ensembl_to_NCBI_df})
            pedigree=$(realpath {input.pedigree})
            mavedb_tsv=$(realpath {params.mavedb_tsv})
            small_report=$(realpath {input.small_variant_report})
            high_impact_report=$(realpath {input.wgs_high_impact_variant_report})
            denovo_report=$(realpath {input.wgs_denovo_variant_report})
            sv_report=$(realpath {input.SV_report})
            cnv_report=$(realpath {input.CNV_report})
            mkdir -p "$stage/input/wgs-coding" "$stage/input/wgs-high-impact" "$stage/input/wgs-denovo" "$stage/reports"
            ln -sf "$small_report" "$stage/input/wgs-coding/{wildcards.family}.wgs.coding.slivar.hg38.csv"
            ln -sf "$high_impact_report" "$stage/input/wgs-high-impact/{wildcards.family}.wgs.high.impact.slivar.hg38.csv"
            ln -sf "$denovo_report" "$stage/input/wgs-denovo/{wildcards.family}.wgs.denovo.slivar.hg38.csv"
            {params.stage_hpo_panel_reports}
            cd "$stage"
            python3 {params.crg2_pacbio}/scripts/annotate_compound_hets.py \
            --seq_type {params.seq_type} \
            --high_med "$high_med" \
            --low "$low" \
            --sv "$sv_report" \
            --cnv "$cnv_report" \
            --ensembl "$ensembl" \
            --ensembl_to_NCBI_df "$ensembl_to_ncbi" \
            --pedigree "$pedigree" \
            {params.hpo_panel_args} \
            --sequence_variant_report_dir input/wgs-coding \
            --wgs_high_impact_variant_report_dir input/wgs-high-impact \
            --wgs_denovo_variant_report_dir input/wgs-denovo \
            --sample_order "$sample_order" \
            --family {wildcards.family} \
            --acmg_sf {params.acmg_sf_flag}
            python3 {params.crg2_pacbio}/scripts/add_mavedb_columns.py \
            --family {wildcards.family} \
            --reports-dir reports \
            --mavedb-tsv "$mavedb_tsv"
            mkdir -p "$root/reports_slivar"
            cp -a reports/. "$root/reports_slivar/") > {log} 2>&1
            """
else:
    rule identify_compound_hets_slivar:
        input:
            high_med_variants="small_variants_slivar/{family}.HIGH-MED.impact.variants.tsv",
            low_variants="small_variants_slivar/{family}.LOW.impact.variants.tsv",
            small_variant_report="small_variants_slivar/coding/{family}/{family}.coding.slivar.hg38.csv",
            wgs_high_impact_variant_report="small_variants_slivar/wgs-high-impact/{family}/{family}.wgs-high-impact.slivar.hg38.csv",
            SV_report="reports/{family}.sv.csv",
            CNV_report="reports/{family}.cnv.csv",
            ensembl=config["annotation"]["general"]["ensembl"],
            ensembl_to_NCBI_df=config["annotation"]["general"]["ensembl_to_NCBI_df"],
            pedigree=config["run"]["ped"],
            sample_order="small_variants/{family}.sample.order.txt",
            **slivar_hpo_panel_inputs,
        output:
            small_variant_report_CH=output_status("reports_slivar/{family}.wgs.coding.CH.hg38.csv"),
            wgs_high_impact_variant_report_CH=output_status("reports_slivar/{family}.wgs.high.impact.CH.hg38.csv"),
            SV_report_CH=output_status("reports_slivar/{family}.sv.CH.hg38.csv"),
            CNV_report_CH=output_status("reports_slivar/{family}.cnv.CH.hg38.csv"),
            compound_het_status="reports_slivar/{family}.compound.het.status.CH.hg38.csv",
            **slivar_hpo_panel_outputs,
        params:
            crg2_pacbio=config["tools"]["crg2_pacbio"],
            seq_type="short",
            hpo_panel_args=get_slivar_hpo_panel_args,
            stage_hpo_panel_reports=stage_slivar_hpo_panel_reports,
            acmg_sf_flag=str(config["run"].get("acmg_sf", "false")).lower(),
            mavedb_tsv=config["annotation"]["general"]["mavedb_tsv"]
        conda:
            "../envs/str_sv.yaml"
        log:
            "logs/compound_hets/{family}.slivar.identify.compound.hets.log",
        shell:
            """
            (set -e
            root=$(pwd)
            stage=$(mktemp -d "$root/.slivar_ch_work.{wildcards.family}.XXXXXX")
            trap 'rm -rf "$stage"' EXIT
            high_med=$(realpath {input.high_med_variants})
            low=$(realpath {input.low_variants})
            sample_order=$(realpath {input.sample_order})
            ensembl=$(realpath {input.ensembl})
            ensembl_to_ncbi=$(realpath {input.ensembl_to_NCBI_df})
            pedigree=$(realpath {input.pedigree})
            mavedb_tsv=$(realpath {params.mavedb_tsv})
            small_report=$(realpath {input.small_variant_report})
            high_impact_report=$(realpath {input.wgs_high_impact_variant_report})
            sv_report=$(realpath {input.SV_report})
            cnv_report=$(realpath {input.CNV_report})
            mkdir -p "$stage/input/wgs-coding" "$stage/input/wgs-high-impact" "$stage/reports"
            ln -sf "$small_report" "$stage/input/wgs-coding/{wildcards.family}.wgs.coding.slivar.hg38.csv"
            ln -sf "$high_impact_report" "$stage/input/wgs-high-impact/{wildcards.family}.wgs.high.impact.slivar.hg38.csv"
            {params.stage_hpo_panel_reports}
            cd "$stage"
            python3 {params.crg2_pacbio}/scripts/annotate_compound_hets.py \
            --seq_type {params.seq_type} \
            --high_med "$high_med" \
            --low "$low" \
            --sv "$sv_report" \
            --cnv "$cnv_report" \
            --ensembl "$ensembl" \
            --ensembl_to_NCBI_df "$ensembl_to_ncbi" \
            --pedigree "$pedigree" \
            {params.hpo_panel_args} \
            --sequence_variant_report_dir input/wgs-coding \
            --wgs_high_impact_variant_report_dir input/wgs-high-impact \
            --sample_order "$sample_order" \
            --family {wildcards.family} \
            --acmg_sf {params.acmg_sf_flag}
            python3 {params.crg2_pacbio}/scripts/add_mavedb_columns.py \
            --family {wildcards.family} \
            --reports-dir reports \
            --mavedb-tsv "$mavedb_tsv"
            mkdir -p "$root/reports_slivar"
            cp -a reports/. "$root/reports_slivar/") > {log} 2>&1
            """
