import pandas as pd
import logging
import os

def log_message(*message):
    if message:
        for i in message:
            logging.info(i)
            print(i)

def get_value(row, col, default="."):
    if col in row.index and pd.notna(row[col]) and str(row[col]).strip() != "":
        return row[col]
    return default

def clean_sample_name(sample_name):
    return str(sample_name).replace("-", "_")

def collapse_sample_zygosity_genotype_values(row, df, value_type):
    sample_values = []

    for col in df.columns:
        if value_type == "zygosity":
            if col.startswith("Zygosity."):
                sample = col.replace("Zygosity.", "")
            elif col.endswith("_zyg"):
                sample = col[:-4]
            else:
                continue

        elif value_type == "genotype":
            if col.endswith("_GT"):
                sample = col[:-3]
            else:
                continue

        else:
            continue

        value = get_value(row, col)
        sample = clean_sample_name(sample)
        sample_values.append(f"{sample}={value}")

    if len(sample_values) == 0:
        return "."

    return ";".join(sample_values)


def make_variant_key(row, input_report_type):
    if input_report_type in ["wgs.coding.CH", "wgs.high.impact.CH", "wgs.denovo.CH"]:
        return f"{get_value(row, 'Position')}:{get_value(row, 'Ref')}:{get_value(row, 'Alt')}"

    if input_report_type in ["sv.CH", "cnv.CH"]:
        chrom = str(get_value(row, "CHROM"))
        if chrom != "." and not chrom.startswith("chr"):
            chrom = f"chr{chrom}"

        return f"{chrom}:{get_value(row, 'POS')}:{get_value(row, 'END')}:{get_value(row, 'SVTYPE')}"

    return "."

def make_acmg_sf_report_rows(df, family, input_report_type, acmg_col):
    acmg_matches = df[df[acmg_col] != "."].copy()
    report_rows = []

    for _, row in acmg_matches.iterrows():
        if input_report_type == "wgs.coding.CH":
            position = get_value(row, "Position")
            gene = get_value(row, "Gene")
            consequence = get_value(row, "Variation")
            ref = get_value(row, "Ref")
            alt = get_value(row, "Alt")
            end = "."
            svtype = "."
            clinvar = get_value(row, "Clinvar")
            gnomad_af = get_value(row, "Gnomad_af")
            ucsc_link = get_value(row, "UCSC_Link")

        elif input_report_type in ["sv.CH", "cnv.CH"]:
            chrom = str(get_value(row, "CHROM"))
            if chrom != "." and not chrom.startswith("chr"):
                chrom = f"chr{chrom}"

            position = f"{chrom}:{get_value(row, 'POS')}"
            gene = get_value(row, "GENE_NAME")
            consequence = get_value(row, "VARIANT")
            ref = "."
            alt = "."
            end = get_value(row, "END")
            svtype = get_value(row, "SVTYPE")
            clinvar = "."
            gnomad_af = get_value(row, "gnomad_GRPMAX_AF")
            ucsc_link = get_value(row, "UCSC_link")

        else:
            continue

        report_rows.append({
            "POSITION": position,
            "END": end,
            "REF": ref,
            "ALT": alt,
            "SVTYPE": svtype,
            "GENE": gene,
            "ACMG_SF_GENE": get_value(row, acmg_col),
            "CONSEQUENCE": consequence,
            "FAMILY": family,
            "SAMPLE_ZYGOSITIES": collapse_sample_zygosity_genotype_values(row, acmg_matches, "zygosity"),
            "SAMPLE_GENOTYPES": collapse_sample_zygosity_genotype_values(row, acmg_matches, "genotype"),
            "CLINVAR": clinvar,
            "GNOMAD_AF": gnomad_af,
            "UCSC_LINK": ucsc_link,
            "IN_HIGH_IMPACT_REPORT": False,
            "IN_DENOVO_REPORT": False,
            "VARIANT_REPORTED_IN": input_report_type,
            "VARIANT_KEY": make_variant_key(row, input_report_type),
        })

    return pd.DataFrame(report_rows)

def update_flags(acmg_sf_report, df, input_report_type, acmg_col):
    if input_report_type not in ["wgs.high.impact.CH", "wgs.denovo.CH"]:
        return acmg_sf_report

    if len(acmg_sf_report) == 0:
        return acmg_sf_report
    
    acmg_matches = df[df[acmg_col] != "."].copy()

    if len(acmg_matches) == 0:
        log_message(f"No ACMG SF variants in {input_report_type}; no flags updated")
        return acmg_sf_report

    flag_variant_keys = set(
        acmg_matches.apply(lambda row: make_variant_key(row, input_report_type), axis=1)
    )

    if input_report_type == "wgs.high.impact.CH":
        acmg_sf_report["IN_HIGH_IMPACT_REPORT"] = (
            acmg_sf_report["IN_HIGH_IMPACT_REPORT"].astype(bool)
            | acmg_sf_report["VARIANT_KEY"].isin(flag_variant_keys)
        )

    if input_report_type == "wgs.denovo.CH":
        acmg_sf_report["IN_DENOVO_REPORT"] = (
            acmg_sf_report["IN_DENOVO_REPORT"].astype(bool)
            | acmg_sf_report["VARIANT_KEY"].isin(flag_variant_keys)
        )

    return acmg_sf_report


def infer_input_report_type(path):
    basename = os.path.basename(str(path))

    for report_type in [
        "wgs.coding.CH",
        "wgs.high.impact.CH",
        "wgs.denovo.CH",
        "sv.CH",
        "cnv.CH",
    ]:
        if f".{report_type}.SF." in basename:
            return report_type

    raise ValueError(f"Could not infer input report type from path: {path}")


def main(family, input_reports, output_csv, acmg_sf_version):
    logfile = f"logs/report/acmg_sf/{family}.acmg_sf_report.log"
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.DEBUG,
        format="%(asctime)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M",
    )

    acmg_col = f"ACMG_SF_v{acmg_sf_version}"
    report_parts = []
    flag_inputs = []

    for input_report in input_reports:
        input_report_type = infer_input_report_type(input_report)
        df = pd.read_csv(input_report)

        if acmg_col not in df.columns:
            raise ValueError(f"{input_report} does not contain expected column {acmg_col}")

        if input_report_type in ["wgs.coding.CH", "sv.CH", "cnv.CH"]:
            rows = make_acmg_sf_report_rows(
                df=df,
                family=family,
                input_report_type=input_report_type,
                acmg_col=acmg_col,
            )
            if len(rows) > 0:
                report_parts.append(rows)

        elif input_report_type in ["wgs.high.impact.CH", "wgs.denovo.CH"]:
            flag_inputs.append((df, input_report_type))

    if len(report_parts) > 0:
        acmg_report = pd.concat(report_parts, ignore_index=True)
    else:
        acmg_report = pd.DataFrame(columns=[
            "POSITION",
            "END",
            "REF",
            "ALT",
            "SVTYPE",
            "GENE",
            "ACMG_SF_GENE",
            "CONSEQUENCE",
            "FAMILY",
            "SAMPLE_ZYGOSITIES",
            "SAMPLE_GENOTYPES",
            "CLINVAR",
            "GNOMAD_AF",
            "UCSC_LINK",
            "IN_HIGH_IMPACT_REPORT",
            "IN_DENOVO_REPORT",
            "VARIANT_REPORTED_IN",
            "VARIANT_KEY",
        ])

    for df, input_report_type in flag_inputs:
        acmg_report = update_flags(
            acmg_sf_report=acmg_report,
            df=df,
            input_report_type=input_report_type,
            acmg_col=acmg_col,
        )

    acmg_report.to_csv(output_csv, index=False)
    log_message(f"{output_csv} created with {len(acmg_report)} rows")


if __name__ == "__main__":
    family = snakemake.wildcards.family
    input_reports = list(snakemake.input.reports)
    output_csv = snakemake.output.report
    acmg_sf_version = snakemake.params.acmg_sf_version

    main(family, input_reports, output_csv, acmg_sf_version)