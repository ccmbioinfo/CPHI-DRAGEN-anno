import numpy as np
import pandas as pd

samples_tsv = snakemake.input.samples_tsv
dragen_qc_summary = snakemake.input.dragen_qc_summary
selfsm_files = snakemake.input.selfsm
sex_check = snakemake.input.sex_check
ped_check = snakemake.input.ped_check

out_tsv = snakemake.output.tsv

samples = list(snakemake.params.samples)
min_cov = float(snakemake.params.min_mean_coverage)
max_freemix = float(snakemake.params.max_freemix)
unrelated_max = float(snakemake.params.unrelated_max_rel)
firstdeg_min = float(snakemake.params.firstdeg_min_rel)
firstdeg_max = float(snakemake.params.firstdeg_max_rel)


def to_bool_str(val):
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return "NA"
    return "True" if val else "False"


# 1. Mean coverage >= min_cov (DRAGEN average autosomal coverage over genome)
cov_summary = pd.read_csv(dragen_qc_summary, sep="\t", dtype=str)
cov_map = {}
for _, row in cov_summary.iterrows():
    try:
        cov_map[row["Sample"]] = float(row["Avg_autosomal_cov"])
    except (ValueError, TypeError, KeyError):
        cov_map[row["Sample"]] = np.nan

coverage_pass = {}
for sample in samples:
    cov = cov_map.get(sample, np.nan)
    coverage_pass[sample] = np.nan if pd.isna(cov) else (cov >= min_cov)


# 2. Contamination (VerifyBamID FREEMIX) < max_freemix (FREEMIX is a fraction 0-1)
freemix_pass = {sample: np.nan for sample in samples}
for sample, selfsm in zip(samples, selfsm_files):
    df = pd.read_csv(selfsm, sep="\t")
    freemix_cols = [c for c in df.columns if c.strip().upper() == "FREEMIX"]
    if df.empty or not freemix_cols:
        continue
    try:
        freemix = float(df.iloc[0][freemix_cols[0]])
    except (ValueError, TypeError):
        continue
    freemix_pass[sample] = freemix < max_freemix


# 3. Derived sex matches pedigree sex (peddy sex_check.csv error == False)
sex_pass = {sample: np.nan for sample in samples}
sex_df = pd.read_csv(sex_check)
if {"sample_id", "error"}.issubset(sex_df.columns):
    for _, row in sex_df.iterrows():
        sid = str(row["sample_id"])
        if sid not in sex_pass:
            continue
        err = str(row["error"]).strip().lower()
        if err in ("true", "false"):
            sex_pass[sid] = err == "false"


# 4. Relatedness consistent with pedigree (peddy ped_check.csv, pairwise -> per-sample)
# See https://help.connected.illumina.com/emedgene/emedgene-analyze-manual/reviewing_a_case/lab_tab/pedigree-section
# Unrelated individuals i.e parents (if not consanguineous) should have a relatedness percentage < 0.2. 
# Parent-child and sibling-sibling should have relatedness between 40 - 60.
# A sample passes only if every pair it participates in is consistent.
rel_checks = {sample: [] for sample in samples}
ped_df = pd.read_csv(ped_check)
if {"sample_a", "sample_b", "rel", "pedigree_relatedness"}.issubset(ped_df.columns):
    for _, row in ped_df.iterrows():
        sample_a = str(row["sample_a"])
        sample_b = str(row["sample_b"])
        try:
            rel = float(row["rel"])
            expected = float(row["pedigree_relatedness"])
        except (ValueError, TypeError):
            continue

        if expected == 0:
            consistent = rel < unrelated_max
        elif expected == 0.5:
            consistent = firstdeg_min <= rel <= firstdeg_max
        else:
            consistent = abs(rel - expected) <= 0.1

        for sample in (sample_a, sample_b):
            if sample in rel_checks:
                rel_checks[sample].append(consistent)

relatedness_pass = {}
for sample in samples:
    checks = rel_checks[sample]
    relatedness_pass[sample] = np.nan if not checks else all(checks)


checks = [
    ("Mean coverage >= {:g}X".format(min_cov), coverage_pass),
    ("Contamination (VerifyBamID FREEMIX) < {:g}%".format(max_freemix * 100), freemix_pass),
    ("Correct sex", sex_pass),
    ("Relatedness consistent with pedigree", relatedness_pass),
]

out_df = pd.DataFrame({"Sample": samples})
for name, result in checks:
    out_df[name] = [to_bool_str(result.get(sample, np.nan)) for sample in samples]

out_df.to_csv(out_tsv, sep="\t", index=False)
