import argparse
from datetime import date
import pandas as pd
import numpy as np


def pivot_repeat_df(repeat_df):
    # convert repeat dataframe long to wide format
    repeat_pivot = repeat_df.pivot(columns=["SAMPLE"], values=["GT", "SO", "motif_count", "REPCI", "ADSP", "ADFL", "ADIR", "coverage"])
    # collapse multi-level index 
    repeat_pivot = repeat_pivot.reset_index()
    cols = repeat_pivot.columns
    new_cols = []
    for col in cols:
        new_col = col[1] + "_" + col[0] if col[1] != "" else col[0]
        new_cols.append(new_col)
    repeat_pivot.columns = new_cols

    return repeat_pivot

def is_disease(motif_count, gene, threshold):
    try:
        motif_count = [int(count) for count in motif_count.split("/")]
    except: # missing genotype
        return "Missing"
    is_disease = False
    if pd.isna(threshold):
        # threshold is NaN
        is_disease = None
    else:
        for count in motif_count:
            if count == ".":
                continue
            elif gene == "VWA1":
                if count != 2: # 2 copies is benign as per STRchive 
                    is_disease = True
            elif count >= int(threshold):
                is_disease = True

    return is_disease


def main(repeat_tsv, disease_thresholds, output_file):
    repeat_df = pd.read_csv(repeat_tsv, sep="\t", names=["SAMPLE", "CHROM", "POS", "VARID", "REF", "RL", "RU", "GT", "SO", "motif_count", "REPCI", "ADSP", "ADFL", "ADIR", "coverage"])
    repeat_df.set_index(["CHROM", "POS", "VARID", "REF", "RL", "RU"], inplace=True)
    samples = repeat_df["SAMPLE"].unique()
    # convert repeat dataframe long to wide format
    repeat_pivot = pivot_repeat_df(repeat_df)
    # remove non-disease loci
    repeat_pivot = repeat_pivot[~repeat_pivot["VARID"].str.contains("chr")].copy()
    # rename some loci to match disease threshold file
    gene_dict = {"HOXA13_1": "HOXA13-I", "HOXA13_2": "HOXA13-II", "HOXA13_3": "HOXA13-III", "ARX_1": "EIEE1_ARX", "ARX_2": "PRTS_ARX","C9ORF72":"C9orf72"}
    repeat_pivot["VARID"] = repeat_pivot["VARID"].replace(gene_dict)
    disease_thresholds = pd.read_csv(disease_thresholds, sep="\t")
    disease_thresholds.columns = disease_thresholds.columns.str.upper()
    repeat_pivot_with_thresholds = repeat_pivot.merge(disease_thresholds, left_on="VARID", right_on="GENE", how="left")
    # now annotate with disease status
    repeat_pivot_with_thresholds.rename(columns={"DISEASE THRESHOLD": "DISEASE_THRESHOLD"}, inplace=True)
    motif_count_cols = repeat_pivot_with_thresholds.columns[repeat_pivot_with_thresholds.columns.str.contains("motif_count")]
    for col in motif_count_cols:
        repeat_pivot_with_thresholds[f"{col}_DISEASE_PREDICTION"] = repeat_pivot_with_thresholds.apply(
            lambda row: is_disease(
                row[col],
                row["GENE"],
                row.DISEASE_THRESHOLD
            ),
            axis=1,
        )
    disease_pred_cols = [col for col in repeat_pivot_with_thresholds.columns if "PREDICTION" in col and "motif_count" in col]
    for col in disease_pred_cols:
        rename_col = col.replace("_motif_count", "")
        repeat_pivot_with_thresholds.rename(columns={col: rename_col}, inplace=True)
    # format for export
    sample_cols = repeat_pivot_with_thresholds.columns[repeat_pivot_with_thresholds.columns.str.contains("EXP")]
    sample_cols = [col for col in sample_cols if "_SO" not in col and "DISEASE_PREDICTION" not in col] 
    disease_pred_cols = [col for col in repeat_pivot_with_thresholds.columns if "DISEASE_PREDICTION" in col]
    report_cols = (
        [
            "CHROM",
            "POS",
            "REF",
            "RL",
            "RU",
            "VARID",
            "DISORDER",
            "DISEASE_THRESHOLD",
        ] + list(disease_pred_cols) + list(sample_cols)
    )

    repeat_pivot_with_thresholds = repeat_pivot_with_thresholds[report_cols].copy()
    repeat_pivot_with_thresholds["DISEASE_THRESHOLD"] = repeat_pivot_with_thresholds["DISEASE_THRESHOLD"].astype(str)
    repeat_pivot_with_thresholds["DISEASE_THRESHOLD"] = repeat_pivot_with_thresholds["DISEASE_THRESHOLD"].str.replace(".0", "")
    repeat_pivot_with_thresholds.rename(columns={"VARID": "GENE", "REF": "REF_REPEAT_COUNT", "RL": "REF_LEN_BP", "RU": "MOTIF"}, inplace=True)
    repeat_pivot_with_thresholds.replace({np.nan: "."}, inplace=True)
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    output_prefix = output_file.replace(".csv", "")
    repeat_pivot_with_thresholds.to_csv(f"{output_file}", index=False)
    repeat_pivot_with_thresholds.to_csv(f"{output_prefix}.{today}.csv", index=False)

















if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generates a report for known pathogenic repeat expansion loci from a DRAGEN genotyped STRs"
    )
    parser.add_argument("--repeat_tsv", type=str, help="Repeat TSV", required=True)
    parser.add_argument("--disease_thresholds", type=str, help="Repeat loci", required=True)
    parser.add_argument("--output_file", type=str, help="Output filename", required=True)

    args = parser.parse_args()
    repeat_tsv = args.repeat_tsv
    disease_thresholds = args.disease_thresholds
    output_file =args.output_file

    main(
        repeat_tsv,
        disease_thresholds,
        output_file
    )