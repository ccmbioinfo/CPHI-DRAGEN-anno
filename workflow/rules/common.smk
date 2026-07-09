import glob
import pandas as pd
import os
from snakemake.utils import validate
from snakemake.utils import min_version
from datetime import date

min_version("9.16.2")


samples = pd.read_table(config["run"]["samples"], dtype=str).set_index("sample", drop=False)
units = pd.read_table(config["run"]["units"], dtype=str).set_index(["family"], drop=False)
family = config["run"]["family"]

def get_small_variant_vcf(wildcards):
    return units.loc[family, "small_variant_vcf"]
    
def get_repeat_dir(wildcards):
    return units.loc[family, "repeat_VCF_dir"]

def get_sv_vcf(wildcards):
    input_vcf = units.loc[family, "SV_vcf"]
    return input_vcf

def get_cram(wildcards):
    return samples.loc[wildcards.sample, "CRAM"]

def get_str_vcf(wildcards):
    return samples.loc[wildcards.sample, "STR"]

def get_metrics(wildcards):
    return samples.loc[wildcards.sample, "metrics"]

def get_dragen_metrics_files(wildcards):
    metrics_files = samples["metrics"].tolist()
    return metrics_files

def get_cnv_vcf(wildcards):
    input_vcf = units.loc[family, "CNV_vcf"]
    return input_vcf

def get_wrapper_path(*dirs):
    return "file:%s" % os.path.join(workflow.basedir, "wrappers", *dirs)
