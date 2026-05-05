import pandas as pd
import logging
import re
import os
from datetime import date

def log_message(*message):
    if message:
        for i in message:
            logging.info(i)
            print(i)

#search report gene string for exact acmg secondary finding gene name macthes 
def find_acmg_sf_gene_matches(report_gene_string, acmg_sf_genes):

    # return empty list for missing, null, or placeholder values
    if pd.isna(report_gene_string) or str(report_gene_string).strip() in ("", "."):
        return []

    report_string = str(report_gene_string).strip()
    matches = []
# (?![A-Za-z0-9]) ensures there is no alphanumeric character on the boundry of gene string matches (ie, not a longer gene name containing the acmg gene within it)
    for gene in acmg_sf_genes:
        pattern = r'(?<![A-Za-z0-9])' + re.escape(gene) + r'(?![A-Za-z0-9])'
        if re.search(pattern, report_string, re.IGNORECASE):
            matches.append(gene)
    return matches

#create acmg report with variants from SNV, SV, and CNV reports
def make_variant_key(row, input_report_type):
    if input_report_type in ["wgs.coding.CH", "wgs.high.impact.CH", "wgs.denovo.CH"]:
        return f"{get_value(row, 'Position')}:{get_value(row, 'Ref')}:{get_value(row, 'Alt')}"
    if input_report_type in ["sv.CH", "cnv.CH"]:
        chrom = str(get_value(row, "CHROM"))

        if chrom != "." and not chrom.startswith("chr"):
            chrom = f"chr{chrom}"

        return f"{chrom}:{get_value(row, 'POS')}:{get_value(row, 'END')}:{get_value(row, 'SVTYPE')}"
    return "."

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

def make_acmg_sf_report_rows(df, family, input_report_type, acmg_col):
    acmg_matches = df[df[acmg_col] != "."].copy()
    report_rows = []
    for _, row in acmg_matches.iterrows():

        # WGS small variant report columns
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

        # SV and CNV report columns
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
            "position": position,
            "end": end,
            "ref": ref,
            "alt": alt,
            "svtype": svtype,
            "gene": gene,
            "acmg_sf_gene": get_value(row, acmg_col),
            "consequence": consequence,
            "family": family,
            "sample_zygosities": collapse_sample_zygosity_genotype_values(row, acmg_matches, "zygosity"),
            "sample_genotypes": collapse_sample_zygosity_genotype_values(row, acmg_matches, "genotype"),
            "clinvar": clinvar,
            "gnomad_af": gnomad_af,
            "ucsc_link": ucsc_link,
            "in_high_impact_report": ".",
            "in_denovo_report": ".",
            "variant_reported_in": input_report_type,
            "variant_key": make_variant_key(row, input_report_type),
        })

    return pd.DataFrame(report_rows)

def update_acmg_sf_report_flags(df, input_report_type, acmg_col, acmg_sf_report_csv):
    if input_report_type not in ["wgs.high.impact.CH", "wgs.denovo.CH"]:
        return

    if not os.path.exists(acmg_sf_report_csv):
        log_message(
            f"{acmg_sf_report_csv} does not exist yet; cannot update {input_report_type} flags"
        )
        return

    acmg_matches = df[df[acmg_col] != "."].copy()

    if len(acmg_matches) == 0:
        log_message(f"No ACMG SF variants in {input_report_type}; no flags updated")
        return

    flag_variant_keys = set(
        acmg_matches.apply(lambda row: make_variant_key(row, input_report_type), axis=1)
    )

    acmg_sf_report = pd.read_csv(acmg_sf_report_csv)

    if "in_high_impact_report" not in acmg_sf_report.columns:
        acmg_sf_report["in_high_impact_report"] = False

    if "in_denovo_report" not in acmg_sf_report.columns:
        acmg_sf_report["in_denovo_report"] = False

    if input_report_type == "wgs.high.impact.CH":
        acmg_sf_report["in_high_impact_report"] = (
            acmg_sf_report["in_high_impact_report"].astype(bool)
            | acmg_sf_report["variant_key"].isin(flag_variant_keys)
        )

    if input_report_type == "wgs.denovo.CH":
        acmg_sf_report["in_denovo_report"] = (
            acmg_sf_report["in_denovo_report"].astype(bool)
            | acmg_sf_report["variant_key"].isin(flag_variant_keys)
        )

    acmg_sf_report.to_csv(acmg_sf_report_csv, index=False)
    log_message(f"{acmg_sf_report_csv} {input_report_type} flags updated")

def main(family, input_report_type, input_csv, output_csv, acmg_sf_tsv, acmg_sf_version, seq_type):
    logfile = f"logs/report/acmg_sf/{family}.{input_report_type}.acmg_sf.log"
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.DEBUG,
        format="%(asctime)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M",
    )
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    suffix = "csv" if seq_type == "long" else "hg38.csv"
    
    acmg_df = pd.read_csv(acmg_sf_tsv, sep="\t")
    acmg_genes = set(acmg_df["Gene"].dropna())
    log_message(f"Loaded {len(acmg_genes)} genes from ACMG SF gene list.")

    df = pd.read_csv(input_csv)
    log_message(f"Loaded {len(df)} rows from {input_csv}")

    gene_col = None
    for col_name in ["Gene", "GENE_NAME", "GENE", "gene"]:
        if col_name in df.columns:
            gene_col = col_name
            break

    if gene_col is None:
        log_message(f"ERROR: No gene column found in {input_report_type} report. Expected header names: Gene, GENE_NAME, GENE, gene.")
        df[f"ACMG_SF_v{acmg_sf_version}"] = "."
        dated_output_csv = f"reports/{family}.{input_report_type}.SF.{today}.{suffix}"
        df.to_csv(dated_output_csv, index=False)
        try:
            if os.path.islink(output_csv) or os.path.exists(output_csv):
                os.remove(output_csv)
            os.symlink(os.path.basename(dated_output_csv), output_csv)
        except Exception as e:
            log_message(f"Could not create symlink {output_csv} -> {dated_output_csv}: {e}")
        return

    acmg_sf_matches = []

    for gene_string in df[gene_col]:
        matches = find_acmg_sf_gene_matches(gene_string, acmg_genes)
        if matches:
            unique_matches = sorted(set(matches))
            acmg_sf_matches.append(";".join(unique_matches))
        else:
            acmg_sf_matches.append(".")

    df[f"ACMG_SF_v{acmg_sf_version}"] = acmg_sf_matches
    
    #create secondary findings variant report
    acmg_col = f"ACMG_SF_v{acmg_sf_version}"
    acmg_sf_report_csv = f"reports/{family}.acmg_sf_report.csv"

    main_acmg_report_types = ["wgs.coding.CH", "sv.CH", "cnv.CH"]
    flag_only_report_types = ["wgs.high.impact.CH", "wgs.denovo.CH"]

    if input_report_type in main_acmg_report_types:
        acmg_sf_report = make_acmg_sf_report_rows(
            df=df,
            family=family,
            input_report_type=input_report_type,
            acmg_col=acmg_col
        )

        if len(acmg_sf_report) > 0:
            if os.path.exists(acmg_sf_report_csv):
                acmg_sf_report.to_csv(
                    acmg_sf_report_csv,
                    mode="a",
                    header=False,
                    index=False
                )
            else:
                acmg_sf_report.to_csv(acmg_sf_report_csv, index=False)

        log_message(f"{acmg_sf_report_csv} updated with {input_report_type} rows")

    elif input_report_type in flag_only_report_types:
        update_acmg_sf_report_flags(
            df=df,
            input_report_type=input_report_type,
            acmg_col=acmg_col,
            acmg_sf_report_csv=acmg_sf_report_csv
        )
    
    num_rows_matching_ACMG_SF_list = (df[f"ACMG_SF_v{acmg_sf_version}"] != ".").sum()
    log_message(f"{num_rows_matching_ACMG_SF_list} variants impacting ACMG SF v{acmg_sf_version} genes")

    dated_output_csv = f"reports/{family}.{input_report_type}.SF.{today}.{suffix}"
    df.to_csv(dated_output_csv, index=False)
    log_message(f"{dated_output_csv} created")
    
    symlink_path = output_csv
    target_path = dated_output_csv
    try:
        if os.path.islink(symlink_path) or os.path.exists(symlink_path):
            os.remove(symlink_path)
        os.symlink(os.path.basename(target_path), symlink_path)
    except Exception as e:
        log_message(f"Could not create symlink {symlink_path} -> {target_path}: {e}")

if __name__ == "__main__":
    family = snakemake.wildcards.family
    input_report_type = snakemake.wildcards.input_report_type
    input_csv = snakemake.input.report
    output_csv = snakemake.output.report
    acmg_tsv = snakemake.input.acmg_sf_list
    acmg_sf_version = snakemake.params.acmg_sf_version
    seq_type = snakemake.params.seq_type
    main(family, input_report_type, input_csv, output_csv, acmg_tsv, acmg_sf_version, seq_type)