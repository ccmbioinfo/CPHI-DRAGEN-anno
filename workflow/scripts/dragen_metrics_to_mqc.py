import pandas as pd


samples_tsv = snakemake.input.samples_tsv

duplication_out = snakemake.output.duplication
mapping_out = snakemake.output.mapping
qc_summary_out = snakemake.output.qc_summary

# Read sample tsv
samples_df = pd.read_csv(samples_tsv, sep="\t", dtype=str)

# Build a combined dataframe
all_rows = []
for _, row in samples_df.iterrows():
    sample = row["sample"]
    metrics_path = f"{row["DRAGEN_results_dir"]}/{sample}.metrics.tsv"

    metrics_df = pd.read_csv(metrics_path, sep="\t", dtype=str)
    if metrics_df.empty:
        raise ValueError(f"No rows found in metrics TSV: {metrics_path}")

    rec = metrics_df.iloc[0].copy()
    rec["Sample"] = sample
    all_rows.append(rec)

all_metrics = pd.DataFrame(all_rows)

dragen_output_schema = snakemake.params.dragen_output_schema

if dragen_output_schema == "modified":
    # Convert DRAGEN metrics columns to the columns expected below.
    new_numeric_cols = [
        "Total input reads",
        "Mapped reads R1",
        "Mapped reads R2",
        "Mapped bases R1",
        "Mapped bases R2",
        "Q30 bases R1",
        "Q30 bases R2",
        "Total bases R1",
        "Total bases R2",
        "Reads with MAPQ [40:inf)",
        "Reads with MAPQ [30:40)",
        "Reads with MAPQ [10:20)",
        "Reads with MAPQ [ 0:10)",
        "Reads with MAPQ NA (Unmapped reads)",
        "PCT of genome with coverage [ 50x:100x)",
        "PCT of genome with coverage [ 20x: 50x)",
        "PCT of genome with coverage [ 15x: 20x)",
        "PCT of genome with coverage [ 10x: 15x)",
        "Reads without mate sequenced",
        "Singleton reads (itself mapped; mate unmapped)",
        "Paired reads mapped to different chromosomes",
        "Paired reads (itself & mate mapped)",
        "Properly paired reads",
    ]

    for col in new_numeric_cols:
        if col in all_metrics.columns:
            all_metrics[col] = pd.to_numeric(all_metrics[col], errors="coerce")

    # Reconstruct summary columns expected by the rest of this script.
    if {"Mapped reads R1", "Mapped reads R2"}.issubset(all_metrics.columns):
        all_metrics["Mapped reads"] = all_metrics["Mapped reads R1"] + all_metrics["Mapped reads R2"]

    if {"Mapped reads", "Total input reads"}.issubset(all_metrics.columns):
        all_metrics["Mapped reads %"] = all_metrics["Mapped reads"] / all_metrics["Total input reads"] * 100.0

    if {"Mapped bases R1", "Mapped bases R2"}.issubset(all_metrics.columns):
        all_metrics["Mapped bases"] = all_metrics["Mapped bases R1"] + all_metrics["Mapped bases R2"]

    if {"Q30 bases R1", "Q30 bases R2", "Total bases R1", "Total bases R2"}.issubset(all_metrics.columns):
        total_q30 = all_metrics["Q30 bases R1"] + all_metrics["Q30 bases R2"]
        total_bases = all_metrics["Total bases R1"] + all_metrics["Total bases R2"]
        all_metrics["Q30 bases %"] = total_q30 / total_bases * 100.0

    mapq_pct_map = {
        "Reads with MAPQ [40:inf)": "Reads with MAPQ [40:inf) %",
        "Reads with MAPQ [30:40)": "Reads with MAPQ [30:40) %",
        "Reads with MAPQ [10:20)": "Reads with MAPQ [10:20) %",
        "Reads with MAPQ [ 0:10)": "Reads with MAPQ [ 0:10) %",
        "Reads with MAPQ NA (Unmapped reads)": "Reads with MAPQ NA (Unmapped reads) %",
    }

    for count_col, pct_col in mapq_pct_map.items():
        if {count_col, "Total input reads"}.issubset(all_metrics.columns):
            all_metrics[pct_col] = all_metrics[count_col] / all_metrics["Total input reads"] * 100.0

    coverage_20_cols = {
        "PCT of genome with coverage [ 50x:100x)",
        "PCT of genome with coverage [ 20x: 50x)",
    }
    if coverage_20_cols.issubset(all_metrics.columns):
        all_metrics["PCT of genome with coverage [  20x: inf)"] = (
            all_metrics["PCT of genome with coverage [ 50x:100x)"] +
            all_metrics["PCT of genome with coverage [ 20x: 50x)"]
        )

    coverage_10_cols = {
        "PCT of genome with coverage [ 50x:100x)",
        "PCT of genome with coverage [ 20x: 50x)",
        "PCT of genome with coverage [ 15x: 20x)",
        "PCT of genome with coverage [ 10x: 15x)",
    }
    if coverage_10_cols.issubset(all_metrics.columns):
        all_metrics["PCT of genome with coverage [  10x: inf)"] = (
            all_metrics["PCT of genome with coverage [ 50x:100x)"] +
            all_metrics["PCT of genome with coverage [ 20x: 50x)"] +
            all_metrics["PCT of genome with coverage [ 15x: 20x)"] +
            all_metrics["PCT of genome with coverage [ 10x: 15x)"]
        )

    if {"Reads without mate sequenced", "Total input reads"}.issubset(all_metrics.columns):
        all_metrics["Reads without mate sequenced %"] = (
            all_metrics["Reads without mate sequenced"] / all_metrics["Total input reads"] * 100.0
        )

    if {"Singleton reads (itself mapped; mate unmapped)", "Total input reads"}.issubset(all_metrics.columns):
        all_metrics["Singleton reads (itself mapped; mate unmapped) %"] = (
            all_metrics["Singleton reads (itself mapped; mate unmapped)"] / all_metrics["Total input reads"] * 100.0
        )

    if {"Paired reads mapped to different chromosomes", "Total input reads"}.issubset(all_metrics.columns):
        all_metrics["Paired reads mapped to different chromosomes %"] = (
            all_metrics["Paired reads mapped to different chromosomes"] / all_metrics["Total input reads"] * 100.0
        )

    if {"Paired reads (itself & mate mapped)", "Properly paired reads"}.issubset(all_metrics.columns):
        all_metrics["Not properly paired reads (discordant)"] = (
            all_metrics["Paired reads (itself & mate mapped)"] - all_metrics["Properly paired reads"]
        )

    if {"Not properly paired reads (discordant)", "Total input reads"}.issubset(all_metrics.columns):
        all_metrics["Not properly paired reads (discordant) %"] = (
            all_metrics["Not properly paired reads (discordant)"] / all_metrics["Total input reads"] * 100.0
        )

# All columns should be numeric but just confirming
numeric_cols = [
    "Total input reads",
    "Paired reads (itself & mate mapped)",
    "Properly paired reads",
    "Properly paired reads %",
    "Number of duplicate marked reads",
    "Number of duplicate marked reads %",
    "Number of unique reads (excl. duplicate marked reads)",
    "Reads without mate sequenced",
    "Reads without mate sequenced %",
    "Unmapped reads",
    "Unmapped reads %",
    "Mapped reads",
    "Mapped reads %",
    "QC-failed reads",
    "QC-failed reads %",
    "Reads mapping to multiple locations",
    "Reads mapping to multiple locations %",
    "Paired reads mapped to different chromosomes",
    "Paired reads mapped to different chromosomes %",
    "Mapped bases",
    "Total trimmed bases",
    "Insert length: mean",
    "Average autosomal coverage over genome",
    "Median autosomal coverage over genome",
    "PCT of genome with coverage [  20x: inf)",
    "PCT of genome with coverage [  10x: inf)",
    "Estimated sample contamination",
    "Q30 bases %",
    "Reads with MAPQ [40:inf) %",
    "Reads with MAPQ [30:40) %",
    "Reads with MAPQ [10:20) %",
    "Reads with MAPQ NA (Unmapped reads) %",
    "Reads with MAPQ [ 0:10) %",
    "Number of unique reads (excl. duplicate marked reads) %",
    "Not properly paired reads (discordant) %",
    "Singleton reads (itself mapped; mate unmapped) %",
    "Mapped reads to ref-external sequences (PAI or NRD) %",
    "Unmapped reads minus ref-external or filtered or excluded %",
    "Not properly paired reads (discordant)",
    "Singleton reads (itself mapped; mate unmapped)",
    "Mapped reads to ref-external sequences (PAI or NRD)",
    "Unmapped reads minus ref-external or filtered or excluded"
]

for col in numeric_cols:
    if col in all_metrics.columns:
        all_metrics[col] = pd.to_numeric(all_metrics[col], errors="coerce")

# Duplication plot TSV
duplication_df = all_metrics[
    [
        "Sample",
        "Number of unique reads (excl. duplicate marked reads)",
        "Number of duplicate marked reads",
    ]
].copy()

duplication_df = duplication_df.rename(
    columns={
        "Number of duplicate marked reads": "Duplicate reads",
        "Number of unique reads (excl. duplicate marked reads)": "Unique reads",
    }
)

duplication_df.to_csv(duplication_out, sep="\t", index=False)

# Mapping plot TSV
for pct_col, out_col in {
    "Reads with MAPQ [40:inf) %": "MAPQ_40_inf_reads",
    "Reads with MAPQ [30:40) %": "MAPQ_30_40_reads",
    "Reads with MAPQ [10:20) %": "MAPQ_10_20_reads",
    "Reads with MAPQ [ 0:10) %": "MAPQ_0_10_reads",
    "Reads with MAPQ NA (Unmapped reads) %": "MAPQ_NA_unmapped_reads",
}.items():

    all_metrics[out_col] = all_metrics["Total input reads"] * all_metrics[pct_col] / 100.0

mapping_df = all_metrics[
    [
        "Sample",
        "MAPQ_40_inf_reads",
        "MAPQ_30_40_reads",
        "MAPQ_10_20_reads",
        "MAPQ_0_10_reads",
        "MAPQ_NA_unmapped_reads",
    ]
].copy()

mapping_df = mapping_df.rename(
    columns={
        "MAPQ_40_inf_reads": "MAPQ >=40 reads",
        "MAPQ_30_40_reads": "MAPQ 30-40 reads",
        "MAPQ_10_20_reads": "MAPQ 10-20 reads",
        "MAPQ_0_10_reads": "MAPQ 0-10 reads",
        "MAPQ_NA_unmapped_reads": "Unmapped reads",
    }
)

mapping_df.to_csv(mapping_out, sep="\t", index=False)

# Summary QC table TSV
if (
    "Total input reads" in all_metrics.columns
    and "Reads with MAPQ [ 0:10) %" in all_metrics.columns
):
    all_metrics["Reads with MAPQ <=10"] = (
        all_metrics["Total input reads"] * all_metrics["Reads with MAPQ [ 0:10) %"] / 100.0
    )
else:
    all_metrics["Reads with MAPQ <=10"] = pd.NA

qc_summary_df = all_metrics[
    [
        "Sample",
        "Total input reads",
        "Mapped reads",
        "Properly paired reads",
        "Number of duplicate marked reads",
        "Reads with MAPQ <=10",
        "Mapped bases",
        "PCT of genome with coverage [  20x: inf)",
        "Estimated sample contamination",
        "Q30 bases %",
        "Insert length: mean",
        "Average autosomal coverage over genome",
        "Median autosomal coverage over genome",
    ]
].copy()

qc_summary_df = qc_summary_df.rename(
    columns={
        "Total input reads": "Total_reads",
        "Mapped reads": "Mapped_reads",
        "Properly paired reads": "Properly_paired_reads",
        "Number of duplicate marked reads": "Duplicate_reads",
        "Reads with MAPQ <=10": "MQ_le_10_reads",
        "Mapped bases": "Mapped_bases",
        "PCT of genome with coverage [  20x: inf)": "Genome_20x_pct",
        "Estimated sample contamination": "Estimated_contamination",
        "Q30 bases %": "Q30_bases_pct",
        "Insert length: mean": "Insert_length_mean",
        "Average autosomal coverage over genome": "Avg_autosomal_cov",
        "Median autosomal coverage over genome": "Median_autosomal_cov",
    }
)

qc_summary_df.to_csv(qc_summary_out, sep="\t", index=False)
