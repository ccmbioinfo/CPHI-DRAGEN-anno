acmg_sf_input_report_type = [
    "wgs.coding.CH",
    "panel.CH",
    "panel-flank.CH",
    "wgs.high.impact.CH",
    "sv.CH",
    "cnv.CH",
]

if len(children) > 0:
    acmg_sf_input_report_type.append("wgs.denovo.CH")


rule add_acmg_sf_column:
    input:
        report="reports/{family}.{input_report_type}.hg38.csv",
        acmg_sf_list=config["annotation"]["general"]["acmg_sf_list"],
    output:
        report="reports/{family}.{input_report_type}.SF.hg38.csv",
    params:
        acmg_sf_version=config["annotation"]["general"]["acmg_sf_version"],
    wildcard_constraints:
        input_report_type="|".join([t.replace(".", "\\.") for t in acmg_sf_input_report_type]),
    log:
        "logs/report/acmg_sf/{family}.{input_report_type}.SF.log",
    conda:
        "../envs/acmg_sf.yaml"
    script:
        "../scripts/add_acmg_sf_column.py"
