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

def main(family, input_report_type, input_csv, output_csv, acmg_sf_tsv, acmg_sf_version):
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
    suffix = "hg38.csv"
    
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
    main(family, input_report_type, input_csv, output_csv, acmg_tsv, acmg_sf_version)