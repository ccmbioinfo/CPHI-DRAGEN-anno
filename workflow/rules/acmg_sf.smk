acmg_sf_input_report_type = [
    "wgs.coding.CH",
    "wgs.high.impact.CH",
    "sv.CH",
    "cnv.CH",
]
if len(children) > 0:
    acmg_sf_input_report_type.append("wgs.denovo.CH")

rule add_acmg_sf_columns:
    input:
        report="reports/{family}.{input_report_type}.hg38.csv",
        acmg_sf_list=config["annotation"]["general"]["acmg_sf_list"],
    output:
        report="reports/{family}.{input_report_type}.SF.hg38.csv",
    params:
        acmg_sf_version=config["annotation"]["general"]["acmg_sf_version"],
        seq_type="short",
    wildcard_constraints:
        input_report_type="|".join([t.replace(".", "\\.") for t in acmg_sf_input_report_type]),
    log:
        "logs/report/acmg_sf/{family}.{input_report_type}.SF.log",
    conda:
        "../envs/acmg_sf.yaml"
    script:
        "../scripts/add_acmg_sf_columns.py"

rule create_acmg_sf_report:
    input:
        reports=lambda wildcards: expand(
            "reports/{family}.{input_report_type}.SF.hg38.csv",
            family=wildcards.family,
            input_report_type=acmg_sf_input_report_type,
        ),
    output:
        report="reports/{family}.ACMG.SF.hg38.csv",
    params:
        acmg_sf_version=config["annotation"]["general"]["acmg_sf_version"],
    log:
        "logs/report/acmg_sf/{family}.acmg_sf_report.log",
    conda:
        "../envs/acmg_sf.yaml"
    script:
        "../scripts/create_acmg_sf_report.py"

rule add_acmg_sf_columns_slivar:
    input:
        report="reports_slivar/{family}.{input_report_type}.hg38.csv",
        acmg_sf_list=config["annotation"]["general"]["acmg_sf_list"],
    output:
        report="reports_slivar/{family}.{input_report_type}.SF.hg38.csv",
    params:
        acmg_sf_version=config["annotation"]["general"]["acmg_sf_version"],
        seq_type="short",
    wildcard_constraints:
        input_report_type="|".join([t.replace(".", "\\.") for t in acmg_sf_input_report_type]),
    log:
        "logs/report/acmg_sf/slivar/{family}.{input_report_type}.SF.log",
    conda:
        "../envs/acmg_sf.yaml"
    script:
        "../scripts/add_acmg_sf_columns.py"

rule create_acmg_sf_report_slivar:
    input:
        reports=lambda wildcards: expand(
            "reports_slivar/{family}.{input_report_type}.SF.hg38.csv",
            family=wildcards.family,
            input_report_type=acmg_sf_input_report_type,
        ),
    output:
        report="reports_slivar/{family}.ACMG.SF.hg38.csv",
    params:
        acmg_sf_version=config["annotation"]["general"]["acmg_sf_version"],
    log:
        "logs/report/acmg_sf/slivar/{family}.acmg_sf_report.log",
    conda:
        "../envs/acmg_sf.yaml"
    script:
        "../scripts/create_acmg_sf_report.py"
