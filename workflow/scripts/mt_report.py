import pandas as pd
import logging
import numpy as np
import io
import gzip
from datetime import date
import os
import re


def log_message(*message):
    """write message to logfile and stdout"""
    if message:
        for i in message:
            logging.info(i)
            print(i)

def concat_df(df1, df2):
    """Concatenate two dataframes along axis 1 (column)"""

    concatenated_df = pd.concat([df1, df2], axis=1)
    log_message("Successfully joined the two dataframes")
    return concatenated_df

def remove_cols(df):
    """Remove unwanted columns from the report dataframe"""

    # List of columns to be removed from the file
    remove_cols = [
        "TIER",
        "REF DEPTH",
        "TOTAL LOCUS DEPTH",
        "VARIANT QUALITY",
        "QUAL",
        "MQM_INFO",
        "MQMR_INFO",
        "QA_INFO",
        "QR_INFO",
        "SAF_INFO",
        "SAR_INFO",
        "SRF_INFO",
        "SRR_INFO",
        "SBR_INFO",
        "SBA_INFO",
        "POS_FILTER",
        "SBR_FILTER",
        "SBA_FILTER",
        "MQMR_FILTER",
        "AQR_FILTER",
        "GT_FORMAT",
        "QR_FORMAT",
        "AQR_FORMAT",
        "QA_FORMAT",
        "AQA_FORMAT",
        "INFO",
        "FORMAT",
    ]

    # remove columns and store the remaining cols in new_df; skip columns that do not exist (eg TIER)
    # implement errors=ignore" rather than deleting the above fields entirely in case some VCFs have them or mity normalise is added to this pipeline
    new_df = df.drop(remove_cols, axis=1, errors ="ignore")
    log_message("Removed unwanted columns from the dataframe")
    return new_df


def split_cols_by_sample(grouped_df):
    samples_info = {}
    for variant, row in grouped_df.iterrows():
        sample = str.split(row.SAMPLE, ";")
        ALT_DEPTH = str.split(row["ALT DEPTH"], ";")
        SAMPLE_DEPTH = str.split(row["TOTAL SAMPLE DEPTH"], ";")
        VARIANT_HETEROPLASMY = str.split(row["VARIANT HETEROPLASMY"], ";")

        x = {}
        for i in range(0, len(sample)):
            x[f"{sample[i]}.VARIANT HETEROPLASMY"] = VARIANT_HETEROPLASMY[i]
            x[f"{sample[i]}.ALT DEPTH"] = ALT_DEPTH[i]
            x[f"{sample[i]}.TOTAL SAMPLE DEPTH"] = SAMPLE_DEPTH[i]
        samples_info[variant] = x
    df = pd.DataFrame(samples_info).transpose().sort_index(axis=1).fillna("-")
    df["HGVS"] = df.index
    df["SAMPLE"] = grouped_df["SAMPLE"]

    return df


def sort_by_sample(df):
    subset_df = df[
        [
            "HGVS",
            "SAMPLE",
            "ALT DEPTH",
            "TOTAL SAMPLE DEPTH",
            "VARIANT HETEROPLASMY",
        ]
    ]
    grouped_df = subset_df.groupby("HGVS").agg(
        {
            "SAMPLE": lambda x: ";".join(str(i) for i in x),
            "ALT DEPTH": lambda x: ";".join(str(i) for i in x),
            "TOTAL SAMPLE DEPTH": lambda x: ";".join(str(i) for i in x),
            "VARIANT HETEROPLASMY": lambda x: ";".join(str(i) for i in x),
        }
    )

    df = df.drop(
        [
            "SAMPLE",
            "ALT DEPTH",
            "TOTAL SAMPLE DEPTH",
            "VARIANT HETEROPLASMY",
        ],
        1,
    )

    df2 = split_cols_by_sample(grouped_df)

    final = pd.merge(df, df2, on="HGVS", how="outer")
    log_message("Report sorted by samples")

    return final.drop_duplicates(ignore_index=True)

  
def get_vcf_info(vcf,report,samples):
#loop over each sample and create sample depth, vaf, and alt depth columns
    for sample in samples:
        sample_depths = []
        vafs = []
        alt_depths = []

        for row in report.iterrows():
            pos=row[1]["POS"]
            ref=row[1]["REF"]
            alt=row[1]["ALT"]
           
            # if pos, ref and alt match with the respective columns in the VCF then get the sample info for that sample
            variant_match = vcf[(vcf["POS"] == pos) & (vcf["REF"] == ref) & (vcf["ALT"] == alt)]
           
            # if no matching VCF row is found append "." for the missing value
            if variant_match.empty:
                sample_depths.append(".")
                vafs.append(".")
                alt_depths.append(".")
                continue
            
            match_row = variant_match.iloc[0]

            # parse the format field for info field matches rather than hardcoding column positions
            # note that the FORMAT fields are record-specific and can change slightly from row to row (eg. PS) 
            format_fields = str(match_row["FORMAT"]).split(":")
            sample_fields = str(match_row[sample]).split(":")
            sample_map = dict(zip(format_fields, sample_fields))

            depth = sample_map.get("DP", ".")
            #uses the original AF field, not the artificial VAF field 
            vaf = sample_map.get("AF", ".")
            #AD is formatted as a tuple; [ref depth, alt depth]
            ad = sample_map.get("AD", ".")
            
            if ad not in [".", ""] and "," in ad:
                alt_depth = ad.split(",")[1]
            else:
                alt_depth = "."

            sample_depths.append(depth)
            vafs.append(vaf)
            alt_depths.append(alt_depth)
#add these info fields to the dataframe for each row in the sample, with blanks as "."
        report[f"{sample}.TOTAL SAMPLE DEPTH"]=sample_depths
        report[f"{sample}.VARIANT HETEROPLASMY"]=vafs
        report[f"{sample}.ALT DEPTH"]=alt_depths

    return report


def check_sort(vcf,df):
    sample = df.SAMPLE.unique()
    if len(sample) == 1:
        log_message("Only one sample present in report")
        return df
    else:
        log_message("Multiple samples present in report")
        updated_df = sort_by_sample(df)
        updated_df=get_vcf_info(vcf,updated_df,sample)
        return updated_df

def keep_only_pass(report):
    report=report[report["FILTER"]=="PASS"]
    return report



# create one field containing collapsed values across the MITOMAP groups 
# region groups from MITOMAP: RNA, Coding, ControL
# disease association groups from MITOMAP: Confirmed, Disease

def clean_empty_value(value):
    if pd.isna(value):
        return None

    value = str(value).strip()

    if value in {"", "."}:
        return None

    return value

def return_first_value(*values):
#Returns the first non-empty value from a priority list
    for val in values:
        val = str(val).strip()
        if val not in ("", "."):
            return val
    return ""

def combine_labeled_values(row, column_map):
    parts = []

    for label, column_name in column_map.items():
        if column_name in row.index:
            value = clean_empty_value(row[column_name])
            if value is not None:
                parts.append(f"{label}: {value}")

    return "; ".join(parts)

def add_additional_reported_disease_associations(row):
    confirmed_association = str(row.get("MITOMAP CONFIRMED MUTATIONS ASSOCIATEDDISEASE", "")).strip()

    for col in ["MITOMAP MUTATIONS RNA DISEASE", "MITOMAP DISEASE DISEASE"]:
        reported_association = str(row.get(col, "")).strip()

        if reported_association in ("", "."):
            continue

        # clean hyphens from Disease source 
        if col == "MITOMAP DISEASE DISEASE":
            reported_association = re.sub(r'-/-', ' / ', reported_association)
            reported_association = re.sub(r'\|-', ' / ', reported_association)
            reported_association = re.sub(r'(?<=[a-z])-(?=[a-z])', ' ', reported_association)
            reported_association = re.sub(r'(?<=[a-z])-(?=[A-Z][a-z])', ' ', reported_association)

        # only return if meaningfully different from confirmed
        if confirmed_association and reported_association.lower() == confirmed_association.lower():
            return "."

        return reported_association

    return "."

def create_collapsed_mitomap_columns(df):
    priority_fields = {
        "MITOMAP AA CHANGE": [
            "MITOMAP CONFIRMED MUTATIONS AMINOACIDCHANGE",
            "MITOMAP DISEASE AACHANGE",
            "MITOMAP VARIANTS CODING AMINOACIDCHANGE", #least interpretable format
            "MITOMAP MUTATIONS CODING CONTROL AMINOACIDCHANGE",
        ],
        "MITOMAP SIGNIFICANCE SUMMARY": [
            "MITOMAP CONFIRMED MUTATIONS STATUSMITOMAPCLINGEN", #from confirmed table
            "MITOMAP MUTATIONS RNA STATUS", #from RNA table
            "MITOMAP DISEASE DISEASE STATUS", #from coding and control tables
        ],       
        #may differ from GENE/LOCUS field, eg MT-TER variants map to MT-TL1 in MITOMAP
        "MITOMAP LOCUS": [
            "MITOMAP CONFIRMED MUTATIONS LOCUS",
            "MITOMAP MUTATIONS RNA LOCUS",
        ],
								"MITOMAP HOMOPLASMY": [
            "MITOMAP MUTATIONS RNA HOMOPLASMY",
            "MITOMAP DISEASE HOMOPLASMY",
        ],
        "MITOMAP HETEROPLASMY": [
            "MITOMAP MUTATIONS RNA HETEROPLASMY",
             "MITOMAP DISEASE HETEROPLASMY",
        ],
    }

    for new_column, original_cols in priority_fields.items():
        df[new_column] = df.apply(lambda row, cols=original_cols: return_first_value(*[row.get(c) for c in cols]),axis=1)
    
    df["MITOMAP ADDITIONAL REPORTED DISEASE"] = df.apply(add_additional_reported_disease_associations, axis=1)
        
    collapsed_fields = {
        "MITOMAP # REFERENCES": {
            "RNA": "MITOMAP MUTATIONS RNA REFERENCES",
            "Coding": "MITOMAP VARIANTS CODING CURATEDREFERENCES",
        },
    }

    for new_column, column_list in collapsed_fields.items():
        df[new_column] = df.apply(lambda row, cl=column_list: combine_labeled_values(row, cl), axis=1)

    log_message("Created collapsed MITOMAP columns")
    return df

def remove_raw_collapsed_mitomap_columns(df):
    raw_mitomap_columns = [
        "MITOMAP DISEASE AC",
        "MITOMAP DISEASE AF",
        "MITOMAP DISEASE AACHANGE",
        "MITOMAP DISEASE HOMOPLASMY",
        "MITOMAP DISEASE HETEROPLASMY",
        "MITOMAP DISEASE DISEASE",
        "MITOMAP DISEASE DISEASE STATUS",
        "MITOMAP DISEASE HGFL",
        "MITOMAP CONFIRMED MUTATIONS LOCUS",
        "MITOMAP CONFIRMED MUTATIONS LOCUSTYPE",
        "MITOMAP CONFIRMED MUTATIONS ALLELE",
        "MITOMAP CONFIRMED MUTATIONS AMINOACIDCHANGE",
        "MITOMAP CONFIRMED MUTATIONS STATUS MITOMAP_CLINGEN",
        "MITOMAP CONFIRMED MUTATIONS LASTUPDATE",
        "MITOMAP MUTATIONS RNA LOCUS",
        "MITOMAP MUTATIONS RNA DISEASE",
        "MITOMAP MUTATIONS RNA ALLELE",
        "MITOMAP MUTATIONS RNA HOMOPLASMY",
        "MITOMAP MUTATIONS RNA HETEROPLASMY",
        "MITOMAP MUTATIONS RNA STATUS",
        "MITOMAP MUTATIONS RNA REFERENCES",
        "MITOMAP VARIANTS CODING AMINOACIDCHANGE",
        "MITOMAP VARIANTS CODING CURATEDREFERENCES",
    ]

    df = df.drop(columns=raw_mitomap_columns, errors="ignore")
    log_message("removed raw un-collapsed MITOMAP columns")
    return df

def reorder_cols(df):
    """Reorder columns in the report dataframe"""

    colnames = df.columns

    variant_heteroplasmy = [x for x in colnames if x.endswith("VARIANT HETEROPLASMY")]
    alt_depth = [x for x in colnames if x.endswith("ALT DEPTH")]
    total_sample_depth = [x for x in colnames if x.endswith("TOTAL SAMPLE DEPTH")]

    df=keep_only_pass(df)

    col_list = [
        "CHR",
        "POS",
        "REF",
        "ALT",
        "HGVS",
        "GENE/LOCUS",
        "GENE/LOCUS DESCRIPTION",
        "COHORT COUNT",
        "SAMPLE",
    ]

    col_list2 = [
        "MITOMAP AA CHANGE",
        "clinvar_significance",
        "clinvar_status",
        "gnomAD_AC_hom",
        "gnomAD_AC_het",
        "gnomAD_AF_hom",
        "gnomAD_AF_het",
        "gnomAD_max_hl",
        "MITOMAP LOCUS",
        "MITOMAP CONFIRMED DISEASE",
        "MITOMAP ADDITIONAL REPORTED DISEASE",
        "MITOMAP SIGNIFICANCE SUMMARY",
        "MITOMAP HOMOPLASMY",
        "MITOMAP HETEROPLASMY",
        "MITOMAP REFERENCES",
        "MITOMAP FREQUENCY",
        "MITOMAP DISEASE PUBMED IDS",
        "MITOMAP POLYMORPHISMS AC",
        "MITOMAP POLYMORPHISMS AF",
        "MITOMAP POLYMORPHISMS HGFL",
        "MITOMAP MUTATIONS RNA MITOTIP",
        "MITOTIP SCORE",
        "MITOTIP PERCENTILE",
        "MITOTIP QUARTILE",
        "MITOTIP SCORE INTERPRETATION",
        "GB_COUNT",
        "GB_PERCENTAGE",
        "ANTICODON",
        "MGRB FREQUENCY",
        "MGRB FILTER",
        "MGRB AC",
        "MGRB AN",
        "MITOMAP # REFERENCES",
								"MITOMAP STATUS",
        "PHYLOTREE MUT",
        "PHYLOTREE HAPLOTYPE",
    ]

    #reordering is tolerant of any missing columns (e.g. commented out in report-config.yaml)
    desired_cols = col_list + variant_heteroplasmy + alt_depth + total_sample_depth + col_list2
    final_col_list = [c for c in desired_cols if c in df.columns]
    reordered_df = df[final_col_list]
                  

    # replace '.'/'-' with '0' for some columns
    replace_col_values = [ 
        "MITOMAP POLYMORPHISMS AF",
        "MITOMAP POLYMORPHISMS AC",
        "gnomAD_AC_hom",
        "gnomAD_AC_het",
        "gnomAD_AF_hom",
        "gnomAD_AF_het",
        "gnomAD_max_hl",
    ]

    #replacing is tolerant of any missing columns
    for col in replace_col_values:
        if col in reordered_df.columns:
            reordered_df[col] = reordered_df[col].replace(".", 0)
            
    if "MITOMAP DISEASE PUBMED IDS" in reordered_df.columns:
        reordered_df["MITOMAP DISEASE PUBMED IDS"] = reordered_df["MITOMAP DISEASE PUBMED IDS"].replace({"0": "."})
    
    log_message(
        "Replaced . and - with 0 for frequency columns and rearanged the columns in the dataframe"
    )
    return reordered_df

def read_vcf(vcf):
    with gzip.open(vcf, 'r') as f:
        lines = [l for l in f if not l.startswith(b'##')]
        str_lines=[]
        for l in lines:
            str_lines.append(l.decode())

    vcf_df=pd.read_csv(
        io.StringIO(''.join(str_lines)),
        sep='\t'
    )

    return vcf_df

def main(vcf, report, family):
    logfile = f"logs/report/mitochondrial/{family}.mitochondrial.report.log"
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.INFO,
        format="%(asctime)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M",
    )

    report_df = pd.read_excel(report,engine="openpyxl")
    vcf_df=read_vcf(vcf)
    
   
    final_report = remove_cols(report_df)
    final_report = create_collapsed_mitomap_columns(final_report)
    final_report = remove_raw_collapsed_mitomap_columns(final_report)
    final_report = check_sort(vcf_df,final_report)
    final_report = reorder_cols(final_report)

    today = date.today()
    today = today.strftime("%Y-%m-%d")

    final_report.to_csv(
        f"reports/{family}.mito.{today}.csv", index=False
    )
    # create a symlink instead of a new copy for the Snakemake target
    symlink_path = f"reports/{family}.mito.csv"
    target_path = f"reports/{family}.mito.{today}.csv"
    try:
        if os.path.islink(symlink_path) or os.path.exists(symlink_path):
            os.remove(symlink_path)
        os.symlink(os.path.basename(target_path), symlink_path)
    except Exception as e:
        log_message(f"Could not create symlink {symlink_path} -> {target_path}: {e}")
    
    log_message(
        "Final formatted report containing annotated list of mitochondrial variants created!"
    )


if __name__ == "__main__":
    family = snakemake.wildcards.family
    vcf= snakemake.input.vcf
    report=snakemake.input.report
    main(vcf,report,family)
