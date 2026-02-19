import pandas as pd
import os
from snakemake.utils import validate
from snakemake.utils import min_version
from datetime import date

min_version("9.16.2")


samples = pd.read_table(config["run"]["samples"], dtype=str).set_index("sample", drop=False)
units = pd.read_table(config["run"]["units"], dtype=str).set_index(["family"], drop=False)
family = config["run"]["family"]

def get_repeat_dir(wildcards):
    return units.loc[family, "repeat_VCF_dir"]