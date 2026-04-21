import re
import pandas as pd
import numpy as np
from pysam import VariantFile
import argparse
from annotation.annotate import prepare_OMIM, annotate_OMIM
from collections import defaultdict
from pybedtools import BedTool
from sigfig import round
from datetime import date
import io
import sys

def rename_SV_cols(annotsv_df):
    annotsv_df.rename(
        columns={
            "SV_chrom": "CHROM",
            "SV_start": "POS",
            "SV_end": "END",
            "SV_length": "SVLEN",
            "SV_type": "SVTYPE",
        },
        inplace=True,
    )
    #normalize SVTYPE like <DEL> to DEL
    annotsv_df["SVTYPE"] = annotsv_df["SVTYPE"].astype(str).str.replace("<", "", regex=False).str.replace(">", "", regex=False)
    return annotsv_df


def annotsv_df_to_bed(annotsv_df):
    annotsv_df = annotsv_df[["CHROM", "POS", "END", "SVTYPE", "SVLEN", "ID"]]
    bed = BedTool.from_dataframe(annotsv_df)
    return bed


def filter_length(row):
    length = row.SVLEN
    if pd.isna(length):
        return True
    elif abs(length) >= 50:
        return True
    else:
        return False


def apply_filter_length(annotsv_df):
    annotsv_df_filtered = annotsv_df.apply(lambda row: filter_length(row), axis=1)
    return annotsv_df_filtered


def filter_benign(row):
    svtype = row.SVTYPE
    loss_AF_max = row.B_loss_AFmax
    gain_AF_max = row.B_gain_AFmax
    ins_AF_max = row.B_ins_AFmax
    inv_AF_max = row.B_inv_AFmax
    if svtype == "DEL":
        if loss_AF_max == "nan":
            return True
        else:
            return False
    elif svtype == "DUP":
        if gain_AF_max == "nan":
            return True
        else:
            return False
    elif svtype == "INS":
        if ins_AF_max == "nan":
            return True
        else:
            return False
    elif svtype == "INV":
        if inv_AF_max == "nan":
            return True
        else:
            return False
    # no benign annotations for BNDs
    else:
        return True


def apply_filter_benign(annotsv_df):
    annotsv_df_filtered = annotsv_df.apply(lambda row: filter_benign(row), axis=1)
    return annotsv_df_filtered


def merge_full_split_annos(annotsv_df):
    def join_str(s):
        vals = [str(x) for x in s.tolist() if pd.notna(x)]
        return ";".join(vals)
    annotsv_split = annotsv_df[annotsv_df["Annotation_mode"] == "split"]
    annotsv_full = annotsv_df[annotsv_df["Annotation_mode"] == "full"]
    DDD_cols = ["DDD_status", "DDD_mode", "DDD_consequence", "DDD_disease", "DDD_pmid"]
    GenCC_cols = ["GenCC_disease", "GenCC_moi", "GenCC_classification", "GenCC_pmid"]
    transcript_cols = [
        "Tx",
        "Tx_start",
        "Tx_end",
        "Overlapped_tx_length",
        "Overlapped_CDS_length",
        "Overlapped_CDS_percent",
        "Frameshift",
        "Exon_count",
        "Location",
        "Location2",
        "Dist_nearest_SS",
        "Nearest_SS_type",
        "Intersect_start",
        "Intersect_end",
    ]
    cols_of_interest = DDD_cols + GenCC_cols + transcript_cols
    annotsv_split_agg = (
        annotsv_split.groupby("AnnotSV_ID")
        .agg(
            {
                "GenCC_disease": join_str,
                "GenCC_moi": join_str,
                "GenCC_classification": join_str,
                "GenCC_pmid": join_str,
                "Tx": join_str,
                "Tx_start": join_str,
                "Tx_end": join_str,
                "Overlapped_tx_length": join_str,
                "Overlapped_CDS_length": join_str,
                "Overlapped_CDS_percent": join_str,
                "Frameshift": join_str,
                "Exon_count": join_str,
                "Location": join_str,
                "Location2": join_str,
                "Nearest_SS_type": join_str,
                "Dist_nearest_SS": join_str,
                "Intersect_start": join_str,
                "Intersect_end": join_str,
                "DDD_status": join_str,
                "DDD_mode": join_str,
                "DDD_consequence": join_str,
                "DDD_disease": join_str,
                "DDD_pmid": join_str,
            }
        )
        .reset_index()
    )

    annotsv_full = annotsv_full[
        [col for col in annotsv_full.columns if col not in cols_of_interest]
    ]
    annotsv_merge = annotsv_split_agg.merge(
        annotsv_full, how="outer", on="AnnotSV_ID"
    ).reset_index()
    return annotsv_merge


def add_hpo(hpo, gene):
    try:
        genes = [g for g in re.split('[;&]', gene) if g]
    except AttributeError:
        return "NA"
    terms = []
    for gene in genes:
        # split by - for intergenic variants, which are annotated as <upstream_gene>-<downstream_gene>
        gene = gene.split("-")
        for g in gene:
            try:
                term = str(hpo[hpo["Gene ID"] == g]["Features"].values[0])
                term = term.replace("; ", ";").split(";")
                term = list(set(term))
                for t in term:
                    if t not in terms:
                        terms.append(t)
            except IndexError:
                pass
    if len(terms) == 0:
        return "nan"
    else:
        terms = ",".join(terms)
        return terms


def add_omim(omim_df, gene):
    try:
        genes = [g for g in re.split('[;&]', gene) if g]
    except TypeError:
        print("Omim")
        print(gene)
        return [".","."]
    phenos = []
    inheritance = []
    for gene in genes:
        # split by - for intergenic variants, which are annotated as <upstream_gene>-<downstream_gene>
        gene = gene.split("-")
        for g in gene:
            try:
                phenos.append(
                    str(omim_df[omim_df["gene_id"] == g]["omim_phenotype"].values[0])
                )
                inheritance.append(
                    str(
                        omim_df[omim_df["gene_id"] == g]["omim_inheritance"].values[0]
                    )
                )
            except IndexError:
                pass
    if len(phenos) != 0:
        phenos = ",".join(phenos)
    else:
        phenos = "."
    if len(inheritance) != 0:
        inheritance = ",".join(inheritance)
    else:
        inheritance = "."
    return [phenos, inheritance]

def get_format_alt_value(sample_field, format_field, tag):
    try:
        keys = format_field.split(":")
        vals = sample_field.split(":")
        idx = keys.index(tag)
        pair = vals[idx].split(",")
        return pair [1] if len(pair) > 1 else "."
    except Exception:
        return "."

def get_format_value(sample_field, format_field, tag):
    try:
        keys = format_field.split(":")
        vals = sample_field.split(":")
        return vals[keys.index(tag)]
    except Exception:
        return "."

def get_genotype(sample_GT):
    genotype = sample_GT.split(":")[0]
    if genotype == "0/1" or genotype == "0|1" or genotype == "1|0":
        zyg = "het"
    elif genotype == "1/1":
        zyg =  "hom"
        genotype = "1//1"
    elif genotype == "./.":
        zyg = genotype
    elif genotype == "0/0":
        zyg =  "-"
    else:
        zyg = genotype
    return [zyg, genotype]

def get_cnv_zygosity(sample_GT, chrom):
    genotype = sample_GT.split(":")[0]
    chrom = str(chrom).replace("chr", "").upper()
    if genotype in {"./.", "./."}:
        return "-"
    if genotype in {"./1", "./2", "1/.", "2/."}:
        return "het"
    
    if chrom in {"X", "Y"}:
        if genotype in {"0", "0/.", "./0"}:
            return "-"
        if genotype in {"1", "1/.", "./1"}:
            return "hemi"
    return get_genotype(sample_GT)[0]

# def get_CN(sample_GT):
    # copy number: for CNV calls
#    CN = sample_GT.split(":")[1]
#    return CN

def get_CN(sample_field, format_field):
    # copy number: for CNV calls 
    # find CN in DRAGEN VCF using the format field rather than using GT
    return get_format_value(sample_field, format_field, "CN")

def get_exon_counts(annotsv_df, exon_bed):
    # original functions from SVRecords/SVAnnotator.py
    # this function includes the PB SV ID, because there can be multiple SVs with the same CHROM POS END SVTYPE that have different ALTs
    # these are uniquely identified by the PB SV ID
    exon_counts = defaultdict(
        int
    )  # incrementing dict and setting value in df is faster than incrementing values in the df
    exon_ref = BedTool(exon_bed)
    sample_bedtool = BedTool(
        list(annotsv_df.reset_index()[["CHROM", "POS", "END", "SVTYPE", "ID"]].values)
    )

    for interval in sample_bedtool.intersect(exon_ref, wa=True):
        exon_counts[
            (
                str(interval.chrom),
                str(interval.start),
                str(interval.stop),
                str(interval[3]),
                str(interval[4]),
            )
        ] += 1

    count_df = pd.Series(exon_counts).to_frame()
    count_df = count_df.reset_index()
    count_df.columns = ["CHROM", "POS", "END", "SVTYPE", "ID", "EXONS_SPANNED"]
    for col in ["POS", "END"]:
        count_df[col] = count_df[col].astype(int)
    annotsv_df = pd.merge(
        annotsv_df, count_df, how="left", on=["CHROM", "POS", "END", "SVTYPE", "ID"]
    ).fillna(value={"EXONS_SPANNED": 0})
    return annotsv_df


def annotate_pop_svs(annotsv_df, pop_svs, cols, variant_type):
    # intersect annotsv and population SV bed file
    # separate out INS and BND calls for window-based matching, DEL, DUP, and INV calls for SVLEN-based matching
    annotsv_INS_BND_bed = annotsv_df_to_bed(annotsv_df[(annotsv_df['SVTYPE'] == 'INS') | (annotsv_df['SVTYPE'] == 'BND')])
    annotsv_DEL_DUP_INV_bed = annotsv_df_to_bed(annotsv_df[(annotsv_df['SVTYPE'] == 'DEL') | (annotsv_df['SVTYPE'] == 'DUP') | (annotsv_df['SVTYPE'] == 'INV')])

    pop_svs = pd.read_csv(pop_svs, sep="\t")
    pop_svs["CHROM"] = pop_svs["CHROM"].astype(str).str.replace("^chr", "", regex=True)

    pop_required_cols = ["CHROM", "POS", "END", "SVTYPE", "SVLEN"] + cols
    if set(pop_required_cols).issubset(pop_svs.columns):
        pop_svs = pop_svs[pop_required_cols].copy()

    pop_svs_INS_BND = pop_svs[(pop_svs['SVTYPE'] == 'INS') | (pop_svs['SVTYPE'] == 'BND')]
    pop_svs_INS_BND_bed = BedTool.from_dataframe(pop_svs_INS_BND)
    pop_svs_DEL_DUP_INV = pop_svs[(pop_svs['SVTYPE'] == 'DEL') | (pop_svs['SVTYPE'] == 'DUP') | (pop_svs['SVTYPE'] == 'INV')]
    pop_svs_DEL_DUP_INV_bed = BedTool.from_dataframe(pop_svs_DEL_DUP_INV)

   # look for population DEL, DUP, and INV calls with 80% reciprocal overlap with sample DEL, DUP, and INV calls
    intersect_cols = [
        "CHROM",
        "POS",
        "END",
        "SVTYPE",
        "SVLEN",
        "ID",
        "CHROM_pop",
        "POS_pop",
        "END_pop",
        "SVTYPE_pop",
        "SVLEN_pop"
    ] + cols
    intersect_DEL_DUP_INV = annotsv_DEL_DUP_INV_bed.intersect(
        pop_svs_DEL_DUP_INV_bed, wa=True, wb=True, F=0.8, f=0.8
    ).to_dataframe(names=intersect_cols)


    if variant_type == "SV":
        # look for population INS/BND within 50 bp of sample INS/BND 
        window_INS_BND = annotsv_INS_BND_bed.window(pop_svs_INS_BND_bed, w=50).to_dataframe(names=intersect_cols)
        if len(window_INS_BND) > 0:
            for c in ["SVLEN", "SVLEN_pop"]:
                window_INS_BND[c] = pd.to_numeric(window_INS_BND[c], errors="coerce")
            # only keep matches where the size fraction is greater than 0.8, e.g.  an insertion of 1000bp will not be matched to a 100bp insertion at the same position
            window_INS_BND["SVLEN"] = window_INS_BND["SVLEN"].abs()
            window_INS_BND["size_fraction"] = window_INS_BND.apply(
                lambda x: (min([x["SVLEN"], x["SVLEN_pop"]])) / (max([x["SVLEN"], x["SVLEN_pop"]])),
                axis=1,
            )
            window_INS_BND = window_INS_BND[
                (window_INS_BND["size_fraction"] >= 0.8) | (window_INS_BND["SVTYPE"] == "BND")
            ]
            window_INS_BND = window_INS_BND.drop(columns=["size_fraction"])

            # now concatenate the window and SVLEN-based matches
            intersect = pd.concat([window_INS_BND, intersect_DEL_DUP_INV])
        else: # may be no matches if population database does not include INS/BND calls (e.g. DGV)
            intersect = intersect_DEL_DUP_INV
    else:
        intersect = intersect_DEL_DUP_INV
    
    # popSV and sample SV must be same type
    if intersect.empty or "SVTYPE" not in intersect.columns or "SVTYPE_pop" not in intersect.columns:
        # no overlaps found; add empty annotation cols and return
        # if intersection is empty or missing expected columns, return original DF with placeholder annotation columns
        annotsv_pop_svs = annotsv_df.copy()
        for col in cols:
            if col not in annotsv_pop_svs.columns:
                annotsv_pop_svs[col] = "."
        return annotsv_pop_svs

    intersect = intersect[intersect["SVTYPE"] == intersect["SVTYPE_pop"]]

    # make a column with SV details, e.g DEL:1:25266309-25324509
    pop_name = cols[0].split("_")[0]
    intersect[f"{pop_name}_SV"] = (
        intersect["SVTYPE_pop"].astype(str)
        + ":"
        + intersect["CHROM_pop"].astype(str)
        + ":"
        + intersect["POS_pop"].astype(str)
        + "-"
        + intersect["END_pop"].astype(str)
    )
    cols.append(f"{pop_name}_SV")
    intersect = intersect[["CHROM", "POS", "END", "SVTYPE", "ID"] + cols]
    # round AFs
    try:
        AF_col = [col for col in cols if "AF" in col][0]
        intersect[AF_col] = [
            round(float(af), sigfigs=3) for af in intersect[AF_col].values
        ]
    except IndexError:
        pass
    # group by SV, joining annotation columns
    intersect = intersect.astype(str)
    intersect = (
        intersect.groupby(["CHROM", "POS", "END", "SVTYPE", "ID"])[cols]
        .agg({col: "; ".join for col in cols})
        .reset_index()
    )
    # get max allele frequency
    try:
        AF_col = [col for col in cols if "AF" in col][0]
        intersect[f"{pop_name}_maxAF"] = intersect[AF_col].apply(
            lambda x: max([float(af) for af in x.split("; ")])
        )
        cols.append(f"{pop_name}_maxAF")
        hom_cols = [col for col in cols if "HOM" in col or "nhomalt" in col]
        if len(hom_cols) > 0:
            hom_col = hom_cols[0]
            intersect[f"{hom_col}_max"] = intersect[hom_col].apply(
                lambda x: max([float(hom) for hom in x.split("; ")])
            )
            cols.append(f"{hom_col}_max")

    # get max allele counts for C4R
    except IndexError:
        try:
            count_cols = ["C4R_AC", "seen_in_C4R_count", "C4R_nhomalt"]
            for col in count_cols:
                intersect[f"{col}_max"] = intersect[col].apply(
                    lambda x: max([int(ac) for ac in x.split("; ")])
                )
                cols.append(f"{col}_max")
        # get max allele counts for TG
        except KeyError:
            count_cols = ["TG_AC", "seen_in_TG_count", "TG_nhomalt"]
            for col in count_cols:
                intersect[f"{col}_max"] = intersect[col].apply(
                    lambda x: max([int(ac) for ac in x.split("; ")])
                )
                cols.append(f"{col}_max")

    # merge population AF dataframe with annotSV df
    for col in ["POS", "END"]:
        intersect[col] = intersect[col].astype(int)
    annotsv_pop_svs = pd.merge(
        annotsv_df,
        intersect,
        how="left",
        on=["CHROM", "POS", "END", "SVTYPE", "ID"],
    ).fillna(value={col: 0 for col in cols})
    max_cols = [col for col in cols if col.endswith("_max")]
    for col in max_cols:
        annotsv_pop_svs[col] = annotsv_pop_svs[col].astype(int)
    return annotsv_pop_svs


def annotate_UCSC(chr, pos, end):
    UCSC_base_URL = (
        '=HYPERLINK("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position='
    )
    UCSC_full_URL = f'{UCSC_base_URL}{chr}:{pos}-{end}","UCSC_link")'
    return UCSC_full_URL


def annotate_repeats(annotsv_df, repeats, variant_type):
    #  sample SV must be fully encompassed by repeat
    annotsv_bed = annotsv_df_to_bed(annotsv_df)
    repeats = pd.read_csv(repeats, sep="\t", header=None, names=["CHROM", "POS", "END", "Adotto_tandem_repeat"])
    repeats["CHROM"] = repeats["CHROM"].astype(str).apply(lambda x: x.replace("chr", ""))
    repeats_bed = BedTool.from_dataframe(repeats)
    intersect = annotsv_bed.intersect(
        repeats_bed,
        wa=True,
        wb=True,
        f=1,
    ).to_dataframe()
    if variant_type == "SV":
        intersect.columns = [
            "CHROM",
            "POS",
            "END",
            "SVTYPE",
            "SVLEN",
            "ID",
            "CHROM_repeat",
            "POS_repeat",
            "END_repeat",
            "Adotto_tandem_repeat",
        ]
        intersect = intersect[
            [
                "CHROM",
                "POS",
                "END",
                "SVTYPE",
                "ID",
                "Adotto_tandem_repeat",
            ]
        ]
    
        # merge repeat dataframe with annotSV df
        for col in ["POS", "END"]:
            intersect[col] = intersect[col].astype(int)
        annotsv_repeat_svs = pd.merge(
            annotsv_df,
            intersect,
            how="left",
            on=["CHROM", "POS", "END", "SVTYPE", "ID"],
        ).fillna(value={"Adotto_tandem_repeat": "."})
        return annotsv_repeat_svs
    else: 
        annotsv_df["Adotto_tandem_repeat"] = "."
        return annotsv_df


def vcf_to_df(vcf_path):
    with open(vcf_path, "r") as f:
        lines = [l for l in f if not l.startswith("##")]

    df = pd.read_csv(
        io.StringIO("".join(lines)),
        dtype={
            "#CHROM": str,
            "POS": int,
            "ID": str,
            "REF": str,
            "ALT": str,
            "QUAL": str,
            "FILTER": str,
            "INFO": str,
        },
        sep="\t",
    ).rename(columns={"#CHROM": "CHROM"})

    # Drop DRAGEN reference blocks / non-variant rows (ALT='.')
    df = df[df["ALT"].astype(str) != "."].copy()

    return df

def parse_snpeff(snpeff_df, variant_type):
    svtype_list = []
    end_list = []
    ann_list = []
    pos_list = []
    keep_idx = []
    # parse out svtype, end coordinate, ANN, and position
    for index, row in snpeff_df.iterrows():
        info = row["INFO"].split(";")
        pos = row["POS"]
        # svtype = [i for i in info if "SVTYPE=" in i][0].split("=")[1]
        alt = str(row.get("ALT", ""))

        # Prefer ALT symbolic allele for CNV/SV type normalization (<DEL>, <DUP>, <INV>, ...)
        if alt.startswith("<") and alt.endswith(">"):
            svtype = alt[1:-1].strip()
        else:
            svtype_tag = next((x for x in info if x.startswith("SVTYPE=")), None)
            if svtype_tag:
                svtype = svtype_tag.split("=", 1)[1].strip()
            else:
                continue

        # normalize brackets/upper just in case
        svtype = svtype.replace("<", "").replace(">", "").upper()
        keep_idx.append(index)
        svtype_list.append(svtype)
        try:
            end = [i for i in info if "END=" in i][0].split("=")[1]
            if svtype == "INS":
                raw_pos = int(pos)
                raw_end = int(end)

                ciend_found = False
                ciend_r = 0
                for i in info:
                    if i.startswith("CIEND="):
                        ciend_r = int(i.split("=")[1].split(",")[1])
                        ciend_found = True
                        break

                if ciend_found:
                    end = raw_end + ciend_r + 1
                elif raw_end <= raw_pos:
                    end = raw_end + 1
                else:
                    end = raw_end
                cipos_l = 0
                for i in info:
                    if i.startswith("CIPOS="):
                        cipos_l = int(i.split("=")[1].split(",")[0])
                        break
                if cipos_l != 0:
                    pos = raw_pos + int(cipos_l)
            elif svtype in ("DEL", "DUP", "BND"):
                CIEND = 0
                for i in info:
                    if "CIEND=" in i:
                        CIEND = int(i.split("=")[1].split(",")[1])
                        break
                end = int(end) + int(CIEND)
            elif variant_type == "CNV":
                # END is required; CIEND is optional in DRAGEN CNV records
                end = int(end)
                ciend_tag = next((x for x in info if x.startswith("CIEND=")), None)
                if ciend_tag:
                    ciend_vals = ciend_tag.split("=", 1)[1].split(",")
                    ciend_r = int(ciend_vals[1]) if len(ciend_vals) > 1 else int(ciend_vals[0])
                    end = end + ciend_r

                # CIEND = [i for i in info if "CIEND=" in i][0].split("=")[1].split(",")[1]
                # end = int(end) + int(CIEND)
            elif svtype == "INV":
                raw_pos = int(pos)
                raw_end = int(end)

                cipos_l = 0
                cipos_r = 0
                ciend_r = None

                for i in info:
                    if i.startswith("CIPOS="):
                        cipos_l = int(i.split("=", 1)[1].split(",")[0])
                        cipos_r = int(i.split("=", 1)[1].split(",")[1])
                    elif i.startswith("CIEND="):
                        ciend_r = int(i.split("=", 1)[1].split(",")[1])

                pos = raw_pos + cipos_l
                if ciend_r is not None:
                    end = raw_end + ciend_r
                else:
                    end = raw_end + cipos_r
        except IndexError:
            end = int(pos) + 1
        end_list.append(end)
        try:
            ann = [i for i in info if "ANN=" in i and "SVANN=" not in i][0]
            ann_list.append(ann)
        except IndexError:
            ann_list.append("NA")
        if svtype in ("DEL", "DUP", "BND"):
            CIPOS = 0
            for i in info:
                if "CIPOS=" in i:
                    CIPOS = int(i.split("=")[1].split(",")[0])
                    break
            pos = int(pos) + int(CIPOS)
        elif variant_type == "CNV":
            CIPOS_exist = False
            for i in info:
                if "CIPOS=" in i:
                    CIPOS = i.split("=")[1].split(",")[0]
                    CIPOS_exist = True
                elif CIPOS_exist == False:
                    CIPOS = 0
            pos = int(pos) + int(CIPOS)
        pos_list.append(pos)
    # avoids length mismatch when some rows are skipped
    snpeff_df = snpeff_df.loc[keep_idx].copy().reset_index(drop=True)
    snpeff_df["SVTYPE"] = svtype_list
    snpeff_df["END"] = end_list
    snpeff_df["ANN"] = ann_list
    snpeff_df["POS"] = pos_list
    # from snpeff ANN field, parse out variant type, impact level (e.g. HIGH), gene name, and Ensembl gene ID
    variant_list_all = []
    impact_list_all = []
    gene_list_all = []
    ens_gene_list_all = []
    for line in snpeff_df["ANN"].values:
        variant_list = []
        impact_list = []
        gene_list = []
        ens_gene_list = []
        line = line.split(";")[-1]
        for anno in line.split(","):
            if "LOF=" in anno or "NMD=" in anno:
                pass
            else:
                try:
                    anno = anno.split("|")
                    variant_list.append(anno[1])
                    impact_list.append(anno[2])
                    # for some BNDs, snpeff doesn't annotate against a gene?
                    gene = anno[3]
                    if gene == "":
                        # add transcript instead
                        gene = anno[6]
                    gene_list.append(gene)
                    # same as above
                    ens_gene = anno[4]
                    if ens_gene == "":
                        ens_gene = anno[6]
                    ens_gene_list.append(ens_gene)
                except:
                    for l in [variant_list, impact_list, gene_list, ens_gene_list]:
                        l.append("NA")
        variant_list_all.append(";".join(set(variant_list)))
        impact_list_all.append(";".join(set(impact_list)))
        gene_list_all.append(";".join(set(gene_list)))
        ens_gene_list_all.append(";".join(set(ens_gene_list)))
    snpeff_df["IMPACT"] = impact_list_all
    snpeff_df["GENE_NAME"] = gene_list_all
    snpeff_df["ENSEMBL_GENE"] = ens_gene_list_all
    snpeff_df["VARIANT"] = variant_list_all
    snpeff_df["CHROM"] = [chr.replace("chr", "") for chr in snpeff_df["CHROM"].values]
    snpeff_df = snpeff_df[
        [
            "CHROM",
            "POS",
            "END",
            "ID",
            "SVTYPE",
            "ANN",
            "VARIANT",
            "IMPACT",
            "GENE_NAME",
            "ENSEMBL_GENE",
        ]
    ].astype(str)
    return snpeff_df


def merge_annotsv_snpeff(annotsv_df, snpeff_df):
    merged = annotsv_df.merge(
        snpeff_df, on=["CHROM", "POS", "END", "SVTYPE", "ID"], how="left"
    )
    merged = merged.drop(columns=["ANN"])
    return merged


def calculate_sample_SV_overlap(sample_pos, sample_end, database_pos, database_end):
    sample_len = sample_end - sample_pos
    overlap_start = max(sample_pos, database_pos)
    overlap_end = min(sample_end, database_end)
    overlap = (overlap_end - overlap_start) / float(sample_len) * 100
    overlap_perc = round(float(overlap), sigfigs=3)
    return overlap_perc


def add_clingen(clingen_df, gene, colname):
    try:
        genes = [g for g in re.split('[;&]', gene) if g]
    except TypeError:
        print("ClinGen")
        print(gene)
        return "NA"
    clingen = []
    for gene in genes:
        # split by - for intergenic variants, which are annotated as <upstream_gene>-<downstream_gene>
        gene = gene.split("-")
        for g in gene:
            try:
                clingen.append(
                    str(clingen_df[clingen_df["Gene"] == g][colname].values[0])
                )

            except IndexError:
                pass
    if len(clingen) != 0:
        clingen = ",".join(clingen)
    else:
        clingen = "."
    return clingen

def add_clingen_regions(annotsv_df, clingen_regions_df):
    """
    Annotate SVs against ClinGen regions
    """
    def join_str(s):
        return ";".join([str(x) for x in s.tolist() if pd.notna(x)])
    # parse Genomic Location to extract CHROM, START, END
    clingen_regions_df["CHROM"] = clingen_regions_df["Genomic Location"].str.split(":").str[0].apply(lambda x: x.replace("chr", ""))
    clingen_regions_df["START"] = clingen_regions_df["Genomic Location"].str.split(":").str[1].str.split("-").str[0]
    clingen_regions_df["END"] = clingen_regions_df["Genomic Location"].str.split(":").str[1].str.split("-").str[1]
    clingen_regions_df = clingen_regions_df[clingen_regions_df["CHROM"] != "tbd"]
    
    # prepare clingen regions for BedTools 
    clingen_regions_bed_df = clingen_regions_df[["CHROM", "START", "END", "ISCA Region Name"]].copy()
    clingen_regions_bed_df.rename(columns={"ISCA Region Name": "clingen_region_curation"}, inplace=True)
    
    # convert annotsv_df to BedTool format
    annotsv_bed = annotsv_df_to_bed(annotsv_df)
    
    # convert clingen regions to BedTool format
    clingen_regions_bed = BedTool.from_dataframe(clingen_regions_bed_df)
    
    # perform intersection
    intersect = annotsv_bed.intersect(
        clingen_regions_bed,
        wa=True,
        wb=True,
    ).to_dataframe()
    
    if len(intersect) > 0:
        intersect.columns = [
            "CHROM",
            "POS",
            "END",
            "SVTYPE",
            "SVLEN",
            "ID",
            "CHROM_region",
            "POS_region",
            "END_region",
            "clingen_region_curation",
        ]
        
        # make a column with region details, e.g 1:25266309-25324509
        intersect["clingen_region"] = (
            intersect["CHROM_region"].astype(str)
            + ":"
            + pd.to_numeric(intersect["POS_region"], errors="coerce").fillna(0).astype(int).astype(str)
            + "-"
            + pd.to_numeric(intersect["END_region"], errors="coerce").fillna(0).astype(int).astype(str)
        )

        # calculate percent of sample SV overlapped by region
        intersect["sample_SV_clingen_region_perc_overlap"] = intersect.apply(
            lambda row: calculate_sample_SV_overlap(row["POS"], row["END"], row["POS_region"], row["END_region"]), axis=1
        ).astype(str)
        # group by SV to aggregate multiple overlapping regions
        intersect_agg = intersect.groupby(["CHROM", "POS", "END", "SVTYPE", "ID"]).agg(
            {
                "clingen_region_curation": join_str,
                "clingen_region": join_str,
                "sample_SV_clingen_region_perc_overlap": join_str,
            }
        ).reset_index()
        
        # convert columns to proper types for merging
        for col in ["POS", "END"]:
            intersect_agg[col] = intersect_agg[col].astype(int)
        intersect_agg["CHROM"] = intersect_agg["CHROM"].astype(str)
        
        # merge with annotsv_df
        annotsv_df = pd.merge(
            annotsv_df,
            intersect_agg,
            how="left",
            on=["CHROM", "POS", "END", "SVTYPE", "ID"],
        ).fillna({"clingen_region_curation": ".", "clingen_region": ".", "sample_SV_clingen_region_perc_overlap": "."})
    else:
        # No overlaps found, add empty columns
        annotsv_df["clingen_region_curation"] = "."
        annotsv_df["clingen_region"] = "."
        annotsv_df["sample_SV_clingen_region_perc_overlap"] = "."
    return annotsv_df

def shorten_gene_list(value):
    if pd.isna(value):
        return value
    s = str(value).strip()
    if s in {".", "", "NA", "nan"}:
        return s

    genes = [g.strip() for g in s.split(";") if g.strip()]
    if len(genes) <= 2:
        return s

    return f"{genes[0]};[...truncated...];{genes[-1]}"

def add_BND_structure(svtype, info, alt):
    if svtype == "BND":
        info_extended = ";".join([info, alt])
        return info_extended
    else:
        return info

def annotate_gene_CDS(annotsv_df, ensembl):
    """
    Annotate loci against Ensembl gene CDS
    """
    ensembl = pd.read_csv(ensembl)
    ensembl = ensembl[ensembl["Feature"] == "CDS"][["Chromosome", "Start", "End", "gene_name"]]
    ensembl_bed = BedTool.from_dataframe(ensembl)
    annotsv_bed = annotsv_df_to_bed(annotsv_df)
    intersect_cols = ["CHROM",
            "POS",
            "END",
            "SVTYPE",
            "SVLEN",
            "ID",
            "CHROM_ens",
            "POS_ens",
            "END_ens",
            "GENE_NAME",
            "GENE_ID",
            "GENE_BIOTYPE",
            "FEATURE"]
    intersect = annotsv_bed.intersect(ensembl_bed, wa=True, wb=True).to_dataframe(names=intersect_cols)
    intersect = intersect.drop_duplicates(subset=["GENE_NAME"]).groupby(["CHROM", "POS", "END", "SVTYPE", "ID"]).agg({"GENE_NAME": ";".join}).reset_index()
    intersect.rename(columns={"GENE_NAME": "GENE_NAME_CDS"}, inplace=True)
    annotsv_df = pd.merge(
        annotsv_df,
        intersect,
        how="left",
        on=["CHROM", "POS", "END", "SVTYPE", "ID"],
    ).fillna(value={"GENE_NAME_CDS": "."})
    return annotsv_df


def annotate_breakpoint_gene(annotsv_df, ensembl, location):
    """
    Annotate breakpoint (SV/CNV start or end) against Ensembl genes
    """
    ensembl = pd.read_csv(ensembl)
    ensembl = ensembl[ensembl["Feature"] == "gene"][["Chromosome", "Start", "End", "gene_name"]]
    ensembl_bed = BedTool.from_dataframe(ensembl)
    # create bed files for start or end position
    loc_bed_df = annotsv_df[["CHROM", "POS", "END", "SVTYPE", "ID"]].copy()
    if location == "START ":
        loc_bed_df["START"] = loc_bed_df["POS"] - 1  
        loc_bed_df["END_bed"] = loc_bed_df["POS"] 
    else:
        loc_bed_df["START"] = loc_bed_df["END"] - 1  
        loc_bed_df["END_bed"] = loc_bed_df["END"] 

    loc_bed_df = loc_bed_df[["CHROM", "START", "END_bed", "POS", "END", "SVTYPE", "ID"]] 
    loc_bed = BedTool.from_dataframe(loc_bed_df)
    # intersect with ensembl genes
    intersect_cols = ["CHROM", 
            "START", 
            "END_bed", 
            "POS", 
            "END", 
            "SVTYPE", 
            "ID", 
            "CHROM_ens",
            "POS_ens",
            "END_ens",
            "GENE_NAME",
            "GENE_ID",
            "GENE_BIOTYPE",
            "FEATURE"]
    intersect = loc_bed.intersect(ensembl_bed, wa=True, wb=True).to_dataframe(names=intersect_cols)
    # aggregate genes
    intersect_agg = intersect.groupby(["CHROM", "POS", "END", "SVTYPE", "ID"]).agg({"GENE_NAME": ";".join}).reset_index()
    intersect_agg.rename(columns={"GENE_NAME": f"GENE_NAME_{location}"}, inplace=True)
    # merge with annotsv_df
    for col in ["POS", "END"]:
        intersect_agg[col] = intersect_agg[col].astype(int)
    intersect_agg["CHROM"] = intersect_agg["CHROM"].astype(str)
    annotsv_df = pd.merge(annotsv_df, intersect_agg, how="left", on=["CHROM", "POS", "END", "SVTYPE", "ID"]).fillna(value={f"GENE_NAME_{location}": "."})
    
    return annotsv_df

 
def main(
    df,
    snpeff_df,
    variant_type,
    omim,
    hpo,
    prefix,
    exon_bed,
    gnomad,
    dgv,
    thousandg,
    ensembl,
    repeats,
    clingen_HI,
    clingen_TS,
    clingen_disease,
    clingen_regions,
    samples_tsv
):
    # filter out SVs < 50bp
    df_len = df[apply_filter_length(df)]
    df_len = df_len.astype(str)
    # merge full and split AnnotSV annos
    df_merge = merge_full_split_annos(df_len)
    samples_df = pd.read_csv(samples_tsv, sep="\t", dtype=str)["sample"]
    sample_cols = [s for s in samples_df if s in df.columns]
    if len(sample_cols) == 0:
        print(f"No sample columns from {samples_tsv} found in AnnotSV table.")
        print(f"Expected sample IDs: {samples_df}")
        print(f"Available columns: {list(df.columns)}")
        sys.exit(1)
                     
    # extract genotype and alt allele depth
    for sample in sample_cols:
        if variant_type == "CNV":
            df_merge[f"{sample}_zyg"] = [
                get_cnv_zygosity(row[sample], row["CHROM"]) for index, row in df_merge.iterrows()
            ]
        else:
            df_merge[f"{sample}_zyg"] = [
                get_genotype(row[sample])[0] for index, row in df_merge.iterrows()
            ]
        df_merge[f"{sample}_GT"] = [
            get_genotype(row[sample])[1] for index, row in df_merge.iterrows()
        ]
        if variant_type == "CNV":
            # df_merge[f"{sample}_CN"] = [
            #     get_CN(row[sample]) for index, row in df_merge.iterrows()
            # ]
            df_merge[f"{sample}_CN"] = [
                get_CN(row[sample], row["FORMAT"]) for index, row in df_merge.iterrows()
            ]
        df_merge[f"{sample}_PR_alt"] = [
            get_format_alt_value(row[sample], row["FORMAT"], "PR") for index, row in df_merge.iterrows()
        ]
        df_merge[f"{sample}_SR_alt"] = [
        get_format_alt_value(row[sample], row["FORMAT"], "SR") for index, row in df_merge.iterrows()
        ]
        df_merge[f"{sample}_VF_alt"] = [
        get_format_alt_value(row[sample], row["FORMAT"], "VF") for index, row in df_merge.iterrows()
        ]
        df_merge[f"{sample}_GQ"] = [
            get_format_value(row[sample], row["FORMAT"], "GQ") for index, row in df_merge.iterrows()
        ]
        df_merge[f"{sample}_FS"] = [
            get_format_value(row[sample], row["FORMAT"], "FS") for index, row in df_merge.iterrows()
        ]
        df_merge[f"{sample}_DQ"] = [
            get_format_value(row[sample], row["FORMAT"], "DQ") for index, row in df_merge.iterrows()
        ]
        df_merge[f"{sample}_DN"] = [
            get_format_value(row[sample], row["FORMAT"], "DN") for index, row in df_merge.iterrows()
        ]
    zyg_cols = [col for col in df_merge.columns if "_zyg" in col]
    gt_cols = [col for col in df_merge.columns if "_GT" in col]
    pr_alt_cols = [col for col in df_merge.columns if "_PR_alt" in col]
    sr_alt_cols = [col for col in df_merge.columns if "_SR_alt" in col]
    vf_alt_cols = [col for col in df_merge.columns if "_VF_alt" in col]
    gq_cols = [col for col in df_merge.columns if "_GQ" in col]
    fs_cols = [col for col in df_merge.columns if "_FS" in col]
    dq_cols = [col for col in df_merge.columns if "_DQ" in col]
    dn_cols = [col for col in df_merge.columns if "_DN" in col]
    cn_cols = [col for col in df_merge.columns if "_CN" in col]

    # filter out benign SVs with AF > 1%
    # df_merge_notbenign = df_merge[apply_filter_benign(df_merge)]

    # add snpeff annos
    snpeff_df = parse_snpeff(snpeff_df, variant_type)
    for col in ["POS", "END"]:
        snpeff_df[col] = snpeff_df[col].astype(int)
        df_merge[col] = df_merge[col].astype(int)

    # merge snpeff and annotsv df
    df_merge = merge_annotsv_snpeff(df_merge, snpeff_df)

    # add HPO terms by gene matching
    df_merge["HPO"] = "."
    try:
        df_merge["HPO"] = [
            add_hpo(hpo, gene) for gene in df_merge["ENSEMBL_GENE"].values
        ]
    except:
        print("No HPO terms")

    # add OMIM phenos and inheritance by gene matching
    print("Preparing OMIM data")
    omim_df = prepare_OMIM(f"{omim}/genemap2.txt")

    df_merge.to_csv("df_merge.csv", index=False)
    df_merge["omim_phenotype"] = [
        add_omim(omim_df, gene)[0] for gene in df_merge["ENSEMBL_GENE"].values
    ]
    df_merge["omim_inheritance"] = [
        add_omim(omim_df, gene)[1] for gene in df_merge["ENSEMBL_GENE"].values
    ]

    # add gnomAD SVs
    print("Adding gnomAD SV frequencies")
    gnomad_cols = [
        "gnomad_NAME",
        "gnomad_GRPMAX_AF",
        "gnomad_AC",
        "gnomad_HOM",
    ]
    df_merge = annotate_pop_svs(df_merge, gnomad, gnomad_cols, variant_type)

    # add DGV SVs
    print("Adding DGV SV frequencies")
    dgv_cols = [
        "DGV_AF",
        "DGV_NUM_SAMPLES_TESTED",
    ]
    df_merge = annotate_pop_svs(df_merge, dgv, dgv_cols, variant_type)

    # add 1000G SVs
    print("Adding 1000G SV frequencies")
    thousandg_cols = ["1000G_AF", "1000G_AC", "1000G_nhomalt"]
    df_merge = annotate_pop_svs(df_merge, thousandg, thousandg_cols, variant_type)

    # add exon counts
    df_merge = get_exon_counts(df_merge, exon_bed)

    # add Ensembl CDS annotation
    print("Adding Ensembl CDS annotation")
    df_merge = annotate_gene_CDS(df_merge, ensembl)

    # add start and end gene symbols
    print("Adding start and end gene symbols")
    df_merge = annotate_breakpoint_gene(df_merge, ensembl, "START")
    df_merge = annotate_breakpoint_gene(df_merge, ensembl, "END")

    # add UCSC genome browser URL
    df_merge["UCSC_link"] = [
        annotate_UCSC(chrom, pos, end)
        for chrom, pos, end in zip(df_merge["CHROM"], df_merge["POS"], df_merge["END"])
    ]

    # add PacBio repeats used in repeat expansion finding tool
    df_merge = annotate_repeats(df_merge, repeats, variant_type)

    # add clingen haploinsufficiency and triplosensitivity scores
    df_merge["clingen_HI"] = [
        add_clingen(clingen_HI, gene, "Score") for gene in df_merge["GENE_NAME"].values
    ]
    df_merge["clingen_TS"] = [
        add_clingen(clingen_TS, gene, "Score") for gene in df_merge["GENE_NAME"].values
    ]
    df_merge["clingen_disease"] = [
        add_clingen(clingen_disease, gene, "Disease")
        for gene in df_merge["GENE_NAME"].values
    ]
    df_merge["clingen_classification"] = [
        add_clingen(clingen_disease, gene, "Classification")
        for gene in df_merge["GENE_NAME"].values
    ]
    # add clingen region curations
    df_merge = add_clingen_regions(df_merge, clingen_regions)

    # if SVLEN is greater than 500,000 reduce the size of "GENE_NAME" and "ENSEMBL_GENE" to include only the first and last
    svlen_numeric = pd.to_numeric(df_merge["SVLEN"], errors="coerce")
    long_sv_mask = svlen_numeric.abs() > 500000

    for col in ["GENE_NAME", "ENSEMBL_GENE"]:
        df_merge.loc[long_sv_mask, col] = df_merge.loc[long_sv_mask, col].apply(shorten_gene_list)

    # define columns to be included in report
    hpo_cols = ["HPO"] if isinstance(hpo, pd.DataFrame) else []

    # add BND directionality information to INFO column
    df_merge["INFO_extended"] = df_merge.apply(
        lambda row: add_BND_structure(row["SVTYPE"], row["INFO"], row["ALT"]), axis=1
    )
    df_merge["INFO_extended"] = df_merge["INFO_extended"].apply(
    lambda s: ";".join(
        [x for x in str(s).split(";") if not x.startswith("CONTIG=")]
    )
    )
    df_merge = df_merge.drop(["INFO"], axis=1)
    df_merge = df_merge.rename(columns={"INFO_extended": "INFO"})

    report_columns = (
        [
            "CHROM",
            "POS",
            "END",
            "SVLEN",
            "SVTYPE",
            "ID",
            "INFO",
            "FILTER",
            "GENE_NAME",
            "ENSEMBL_GENE",
            "GENE_NAME_CDS",
            "GENE_NAME_START",
            "GENE_NAME_END",
            "VARIANT",
            "IMPACT",
            "UCSC_link",
            "omim_phenotype",
            "omim_inheritance",
            "clingen_disease",
            "clingen_classification",
            "clingen_HI",
            "clingen_TS",
            "clingen_region_curation",
            "clingen_region",
            "sample_SV_clingen_region_perc_overlap",
        ]
        + hpo_cols
        + zyg_cols
        + gt_cols
        + pr_alt_cols
        + sr_alt_cols
        + vf_alt_cols 
        + gq_cols
        + fs_cols
        + dq_cols
        + dn_cols
        + cn_cols
        + ["Tx", "Frameshift", "EXONS_SPANNED", "Nearest_SS_type", "Dist_nearest_SS"]
        + gnomad_cols
        + dgv_cols
        + thousandg_cols
        + [
            "ExAC_delZ",
            "ExAC_dupZ",
            "ExAC_cnvZ",
            "ExAC_synZ",
            "ExAC_misZ",
            "ExAC_pLI",
            "CytoBand",
            "RE_gene",
            "TAD_coordinate",
            "ENCODE_experiment",
            "Repeat_type_left",
            "Repeat_type_right",
            "SegDup_left",
            "SegDup_right",
            "Adotto_tandem_repeat",
            "ENCODE_blacklist_characteristics_left",
            "ENCODE_blacklist_characteristics_right",
        ]
    )

    df_merge = df_merge[report_columns]
    df_merge.rename(columns={"gnomad_HOM": "gnomad_nhomalt", "gnomad_HOM_max": "gnomad_nhomalt_max"}, inplace=True)
    if variant_type == "CNV":
        # exclude splice site annotations for CNVs
        df_merge = df_merge.drop(columns=["Nearest_SS_type", "Dist_nearest_SS", "ID"] + pr_alt_cols + sr_alt_cols + vf_alt_cols + gq_cols + fs_cols)
    df_merge = df_merge.replace("nan", ".")
    df_merge = df_merge.fillna(".")
    df_merge = df_merge.drop_duplicates()
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    df_merge.to_csv(f"{prefix}.{variant_type.lower()}.csv", index=False)
    df_merge.to_csv(f"{prefix}.{variant_type.lower()}.{today}.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generates a structural variant report using an AnnotSV-annotated pbsv VCF generated from PacBio HiFi WGS"
    )
    parser.add_argument("-annotsv", type=str, help="AnnotSV tsv file", required=True)
    parser.add_argument("-snpeff", type=str, help="Snpeff vcf file", required=True)
    parser.add_argument("-variant_type", type=str, help="Variant type: SV or CNV", required=True)
    parser.add_argument(
        "-hpo",
        help="Tab delimited file containing gene names and HPO terms",
        type=str,
        required=False,
    )
    parser.add_argument(
        "-omim",
        type=str,
        required=True,
        help="Path to directory containing OMIM mim2gene and morbidmap files",
    )
    parser.add_argument(
        "-exon",
        help="exon bed file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-gnomad",
        help="gnomad SVs in bed format",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-dgv",
        help="DGV CNV calls in tsv format",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-thousandg", 
        type=str, 
        required=True, 
        help="1000G SV in bed format"
        )
    parser.add_argument(
        "-ensembl",
        help="Ensembl GTF subset CSV file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-samples",
        help="samples.tsv with a 'sample' col",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-repeats",
        help="Adotto repeats",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-clingen_HI",
        help="clingen HI scores",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-clingen_TS",
        help="clingen TS scores",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-clingen_disease",
        help="clingen disease association",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-clingen_regions",
        help="clingen region curations",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    # pull chrom, pos, end, SVtype, format fields and DDD annotations from AnnotSV text file
    df = pd.read_csv(args.annotsv, sep="\t", low_memory=False)
    # pull gene annotations from SnpEff VCF
    snpeff_df = vcf_to_df(args.snpeff)
    df = rename_SV_cols(df)
    prefix = args.annotsv.replace(".AnnotSV.tsv", "")

    hpo = args.hpo
    if hpo:
        hpo = pd.read_csv(
            args.hpo,
            comment="#",
            skip_blank_lines=True,
            sep="\t",
            encoding="ISO-8859-1",
            engine="python",
        )
        # Phenotips TSV has a space in column name: " Gene symbol"
        hpo.columns = hpo.columns.str.strip()

    clingen_HI = pd.read_csv(
        args.clingen_HI,
        sep="\t",
        comment="#",
        header=None,
        names=["CHROM", "POS", "END", "Gene", "Score"],
    )
    clingen_TS = pd.read_csv(
        args.clingen_TS,
        sep="\t",
        comment="#",
        header=None,
        names=["CHROM", "POS", "END", "Gene", "Score"],
    )
    clingen_disease = pd.read_csv(args.clingen_disease, comment="#")
    clingen_region_cols = ["ISCA ID", "ISCA Region Name", "cytoBand", "Genomic Location", "Haploinsufficiency Score", "Haploinsufficiency Description", "Haploinsufficiency PMID1", "Haploinsufficiency PMID2", "Haploinsufficiency PMID3", "Haploinsufficiency PMID4", "Haploinsufficiency PMID5", "Haploinsufficiency PMID6", "Triplosensitivity Score", "Triplosensitivity Description", "Triplosensitivity PMID1", "Triplosensitivity PMID2", "Triplosensitivity PMID3", "Triplosensitivity PMID4", "Triplosensitivity PMID5", "Triplosensitivity PMID6", "Date Last Evaluated", "Haploinsufficiency Disease ID", "Triplosensitivity Disease ID"]
    clingen_regions = pd.read_csv(args.clingen_regions, comment="#", sep="\t", names=clingen_region_cols)

    main(
        df,
        snpeff_df,
        args.variant_type,
        args.omim,
        hpo,
        prefix,
        args.exon,
        args.gnomad,
        args.dgv,
        args.thousandg,
        args.ensembl,
        args.repeats,
        clingen_HI,
        clingen_TS,
        clingen_disease,
        clingen_regions,
        args.samples
    )
