import numpy as np
import pandas as pd
import pyranges as pr
import re
from typing import Optional

def pivot_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert hits dataframe from long to wide format so each sample is associated with two columns:
    [sample]_allele_len, and [sample]_z_score
    """
    hits = df.copy()
    hits["trid"] = hits["trid"] + "_" + hits["allele_type"]
    hits = hits.drop("allele_type", axis=1)
    hit_pivot = hits.pivot(
        index=["trid", "range", "cutoff", "allele_len_std"],
        columns="sample",
        values=["allele_len", "z_score_len", "z_score_len_rank", "AM", "z_score_AM", "MP", "z_score_MP", "LPS", "z_score_LPS", "LPS_rank"],
    ).reset_index()

    hit_pivot.columns = (
        hit_pivot.columns.get_level_values(1)
        + "_"
        + hit_pivot.columns.get_level_values(0)
    )
    hit_pivot = hit_pivot.rename(
        columns={
            "_trid": "trid",
            "_range": "range",
            "_cutoff": "cutoff",
            "_allele_len_std": "allele_len_std",
            "_AM": "AM",
            "_z_score_AM": "z_score_AM",
            "_MP": "MP",
            "_z_score_MP": "z_score_MP",
            "_z_score_len_rank": "z_score_len_rank",
            "_LPS": "LPS",
            "_z_score_LPS": "z_score_LPS",
            "_z_score_LPS_rank": "LPS_rank"
        }
    )

    return hit_pivot

def hits_to_pr(hits: pd.DataFrame) -> pr.PyRanges:
    """
    Convert hits dataframe to PyRanges object
    """
    hits["Chromosome"] = hits["trid"].map(lambda x: x.split("_")[0].replace("chr", ""))
    hits["Start"] = hits["trid"].map(lambda x: x.split("_")[1])
    hits["End"] = hits["trid"].map(lambda x: x.split("_")[2])
    hit_pr = pr.PyRanges(hits)

    return hit_pr

def annotate_motif(hits: pd.DataFrame, repeat_catalog: str) -> pd.DataFrame:
    """
    Annotate TRID motif. These are stripped from the TRID in the outlier file to save space, especially for large families.
    """
    repeat_catalog = pd.read_csv(repeat_catalog, sep="\t", header=None, names=["chr", "pos", "end", "ID"])
    repeat_catalog["trid_short"] = repeat_catalog["ID"].str.split(";").str[0].str.replace("ID=", "").str.rsplit("_",n=1).str[0]
    repeat_catalog["motif"] = repeat_catalog["ID"].str.split(";").str[1].str.replace("MOTIFS=", "")
    repeat_catalog = repeat_catalog[["trid_short", "motif"]]
    hits["trid_short"] = hits["trid"].str.rsplit("_", n=1).str[0]
    hits = hits.merge(repeat_catalog, on="trid_short", how="left").drop("trid_short", axis=1)

    return hits

def prepare_Ensembl_GTF(gtf_path: str, cols: list) -> pd.DataFrame:
    """
    Convert GTF file into a PyRanges object and clean up columns
    """
    ensembl_gr = pr.read_gtf(gtf_path)[cols]
    introns = ensembl_gr.features.introns(by="gene")
    ensembl_gr = pd.concat([ensembl_gr.df, introns.df], ignore_index=True)
    features = ["exon", "three_prime_utr", "five_prime_utr", "intron", "upstream", "downstream"]
    ensembl_gr = ensembl_gr[ensembl_gr.Feature.isin(features)]

    return ensembl_gr

def prepare_TRGT_loci(loci_path: str) -> pr.PyRanges:
    """
    Convert repeat loci derived from TRGTdb into a PyRanges object
    """
    loci = pd.read_csv(loci_path).rename(
        {"chrom": "Chromosome", "end": "End", "start": "Start"}, axis=1
    )
    # remove chr prefix from chromosomes for compatibility with gene table
    loci["Chromosome"] = loci["Chromosome"].apply(lambda chrom: chrom.strip("chr"))
    loci_gr = pr.PyRanges(loci)

    return loci_gr

def annotate_genes(loci: pr.PyRanges, genes: pr.PyRanges) -> pd.DataFrame:
    """
    Annotate loci against Ensembl genes
    """
    loci_ensembl = loci.join(
        genes, suffix="_ensembl", how="left", apply_strand_suffix=False
    ).drop(["Start_ensembl", "End_ensembl", "Strand"])
    loci_ensembl_df = loci_ensembl.df
    loci_ensembl_df["gene_name"] = np.where(
        loci_ensembl_df["gene_name"].isnull(), "-1", loci_ensembl_df["gene_name"]
    )

    return loci_ensembl_df

def add_constraint(constraint: pd.DataFrame, hits_gene: pd.DataFrame) -> pd.DataFrame:
    """Add gnomAD v4 LOEUF and pLI scores for MANE transcripts"""
    constraint = constraint[constraint["mane_select"] == True]
    hits_gene = hits_gene.merge(
        constraint, left_on="gene_name", right_on="gene", how="left"
    )

    return hits_gene

def group_by_gene(hits_gene: pd.DataFrame) -> pd.DataFrame:
    """
    One hit may be associated with multiple gene features and so multiple rows
    Aggregate by gene and join features to remove duplicate rows
    """
    gene_cols = [
        "gene_name",
        "gene_id",
        "gene_biotype",
        "Feature",
        "lof.oe_ci.upper", 
        "lof.pLI",
        "gene",
        "omim_phenotype", 
        "omim_inheritance",
        "HPO"
    ]

    hits_gene[["lof.oe_ci.upper", "lof.pLI", "gene", "omim_phenotype", "omim_inheritance", "HPO"]] = hits_gene[["lof.oe_ci.upper", "lof.pLI", "gene", "omim_phenotype", "omim_inheritance", "HPO"]].astype(str)
    hits_gene_dedup = hits_gene.groupby(["trid"]).agg(
        {
            "gene_name": ";".join,
            "gene_id": ";".join,
            "gene_biotype": ";".join,
            "Feature": ";".join,
            "lof.oe_ci.upper": lambda x: ";".join(np.unique(x).astype(str)),
            "lof.pLI": lambda x: ";".join(np.unique(x).astype(str)),
            "gene": lambda x: ";".join(np.unique(x).astype(str)),
            "omim_phenotype": lambda x: ";".join(np.unique(x).astype(str)),
            "omim_inheritance":  lambda x: ";".join(np.unique(x).astype(str)),
            "HPO": lambda x: ";".join(np.unique(x).astype(str))
        }
    )
    # merge with original loci table
    hits_gene = hits_gene.drop(gene_cols, axis=1)
    hits_gene_merged = hits_gene.merge(hits_gene_dedup, on=["trid"], how="left")
    hits_gene_merged_dedup = hits_gene_merged.drop_duplicates()

    return hits_gene_merged_dedup


def group_by_segdup(hits_gene: pd.DataFrame) -> pd.DataFrame:
    """
    One hit may be associated with multiple segdup features and so multiple rows
    Aggregate by segdups and join features to remove duplicate rows
    """
    hits_segdup_dedup = hits_gene.groupby(["trid"]).agg(
        {
            "Segdup": ";".join,
        }
    )
    # merge with original loci table
    hits_gene = hits_gene.drop("Segdup", axis=1)
    hits_segdup_merged = hits_gene.merge(hits_segdup_dedup, on=["trid"], how="left")
    hits_segdup_merged_dedup = hits_segdup_merged.drop_duplicates(keep="first")

    return hits_segdup_merged_dedup


def prepare_OMIM(genemap2_path: str) -> pd.DataFrame:
    """
    Merge OMIM mim2gene table with morbidmap table to associate phenotypes to genes.
    See omim.org/downloads for download access to these files.
    """
    genemap2 = pd.read_csv(
        genemap2_path,
        sep="\t",
        comment="#",
        names=[
            "Chrom",
            "Start",
            "End",
            "Cyto",
            "Computed_Cyto",
            "MIM_Number",
            "Gene_Symbols",
            "Gene_Name",
            "Approved_Gene_Symbol",
            "Entrez_Gene_ID",
            "Ensembl_Gene_ID",
            "Comments",
            "Phenotypes",
            "Mouse Gene Symbol/ID",
        ],
    ).set_index("MIM_Number")

    # remove genes that are not associated with any phenotype and rename columns
    genemap2 = genemap2[genemap2["Phenotypes"].notna()].rename(
        {"Ensembl_Gene_ID": "gene_id"}, axis=1
    )
    # extract parsed phenotypes and inheritances
    genemap2_parsed = parse_OMIM_phenotypes(genemap2)
    # subset columns
    genemap2_parsed = genemap2_parsed[["gene_id", "omim_phenotype", "omim_inheritance"]]

    return genemap2_parsed


def parse_OMIM_phenotypes(genemap2: pd.DataFrame) -> pd.DataFrame:
    """
    Parse 'Phenotypes' column in genemap2 to extract two new columns, 'omim_phenotype' and 'omim_inheritance'.
    Adapted from https://github.com/OMIM-org/genemap2-parser/blob/master/parseGeneMap2.py.
    """
    phenotypes = genemap2["Phenotypes"].tolist()
    phenotypes_parsed = []
    inheritances_parsed = []
    inheritance_dict = {
        "Autosomal dominant": "AD",
        "Autosomal recessive": "AR",
        "X-linked": "XL",
        "X-linked dominant": "XLD",
        "X-linked recessive": "XLR",
    }

    # Parse the phenotypes
    for phenotype in phenotypes:
        if not pd.isna(phenotype):
            phenotype_parsed = []
            inheritance_parsed = []
            for phenotype in phenotype.split(";"):
                inheritance = ""
                # Clean the phenotype
                phenotype = phenotype.strip()

                # Long phenotype
                matcher = re.match(r"^(.*),\s(\d{6})\s\((\d)\)(|, (.*))$", phenotype)
                if matcher:
                    # Get the fields
                    phenotype = matcher.group(1)
                    phenotype_parsed.append(phenotype)
                    inheritances = matcher.group(5)

                    # Get the inheritances, may or may not be there
                    if inheritances:
                        inheritance_parsed = []
                        for inheritance in inheritances.split(","):
                            inheritance = inheritance.strip()
                            try:
                                inheritance_abbrev = inheritance_dict[inheritance]
                            except KeyError:
                                inheritance_abbrev = inheritance
                            inheritance_parsed.append(inheritance_abbrev)
                    else:
                        inheritance_parsed.append(inheritance)

                # Short phenotype
                else:
                    # Short phenotype
                    matcher = re.match(r"^(.*)\((\d)\)(|, (.*))$", phenotype)
                    if matcher:
                        # Get the fields
                        phenotype = matcher.group(1)
                        phenotype_parsed.append(phenotype)
                        phenotypeMappingKey = matcher.group(2)
                        inheritances = matcher.group(3)

                        # Get the inheritances, may or may not be there
                        if inheritances:
                            inheritance_parsed = []
                            for inheritance in inheritances.split(","):
                                inheritance = inheritance.strip()
                                try:
                                    inheritance_abbrev = inheritance_dict[inheritance]
                                except KeyError:
                                    inheritance_abbrev = inheritance
                                inheritance_parsed.append(inheritance_abbrev)
                        else:
                            inheritance_parsed.append(inheritance)
            phenotypes_parsed.append((",").join(phenotype_parsed))
            inheritances_parsed.append((",").join(inheritance_parsed))
        else:
            phenotypes_parsed.append("")
            inheritances_parsed.append("")

    genemap2["omim_phenotype"] = phenotypes_parsed
    genemap2["omim_inheritance"] = inheritances_parsed
    genemap2["omim_inheritance"] = genemap2["omim_inheritance"].replace(
        {
            "Autosomal dominant": "AD",
            "Autosomal recessive": "AR",
            "X-linked": "XL",
            "X-linked dominant": "XLD",
            "X-linked recessive": "XLR",
        }
    )

    return genemap2


def annotate_OMIM(loci_ensembl: pd.DataFrame, omim: pd.DataFrame) -> pd.DataFrame:
    """
    Annotate loci against OMIM phenotypes. Loci must already be annotated against Ensembl genes.
    """
    loci_ensembl_omim = loci_ensembl.merge(omim, how="left", on="gene_id")
    loci_ensembl_omim["omim_phenotype"] = loci_ensembl_omim["omim_phenotype"].fillna(
        "."
    )

    return loci_ensembl_omim


def gene_set(genes: str) -> str:
    """
    Reduce semi-colon separated gene identifiers to unique gene identifiers
    """
    genes = ";".join(list(set(genes.split(";"))))

    return genes


def add_hpo(hpo: pd.DataFrame, loci_ensembl: pd.DataFrame) -> list:
    """
    Add gene-based HPO terms
    """
    hpo = hpo[["Gene ID", "Features"]].copy()
    hpo.rename({"Gene ID": "gene_id", "Features": "HPO"}, axis=1, inplace=True)
    loci_ensembl_hpo = loci_ensembl.merge(hpo, how="left", on="gene_id")
    loci_ensembl_hpo["HPO"] = loci_ensembl_hpo["HPO"].fillna("-1")

    return loci_ensembl_hpo


def filter_outliers(row: pd.Series, allele_len_cols: list) -> bool:
    """
    Filter out non-outliers (not expanded in any individual)
    """
    allele_lens = [
        float(row[allele_len]) if not pd.isna(row[allele_len]) else 0
        for allele_len in allele_len_cols
    ]
    cutoff = row["cutoff"]
    try:
        if max(allele_lens) <= round(cutoff):
            return False
        else:
            return True
    except:
        return True # missing in control database


def num_expanded(row: pd.Series, allele_len_cols: list, z_score_cols: list) -> bool:
    """
    Sum number of individuals in whom repeat is expanded (i.e. allele length > cutoff AND Z-score > 3)
    """
    allele_lens = [row[allele_len] for allele_len in allele_len_cols]
    z_scores = [row[z_score] for z_score in z_score_cols]
    cutoff = row["cutoff"] if not pd.isna(row["cutoff"]) else 0
    greater_than = sum([allele_len >= cutoff and z_score >= 3 for allele_len, z_score in zip(allele_lens, z_scores)])

    return greater_than


def read_csv_raw_df(sample_path: str) -> pd.DataFrame:
    """
    Read TRGT-denovo CSV
    Code source: https://github.com/PacificBiosciences/trgt-denovo/blob/main/scripts/python/trio_analysis.ipynb
    """
    dtypes = {
        "sample_id": "object",
        "trid": "object",
        "genotype": "int64",
        "denovo_coverage": "int64",
        "allele_coverage": "int64",
        "allele_ratio": "float64",
        "child_coverage": "int64",
        "child_ratio": "float64",
        "mean_diff_father": "float64",
        "mean_diff_mother": "float64",
        "father_dropout_prob": "float64",
        "mother_dropout_prob": "float64",
        "allele_origin": "object",
        "denovo_status": "object",
        "per_allele_reads_father": "object",
        "per_allele_reads_mother": "object",
        "per_allele_reads_child": "object",
        "father_dropout": "object",
        "mother_dropout": "object",
        "child_dropout": "object",
        "index": "int64",
        "father_MC": "object",
        "mother_MC": "object",
        "child_MC": "object",
        "father_AL": "object",
        "mother_AL": "object",
        "child_AL": "object",
        "father_overlap_coverage": "object",
        "mother_overlap_coverage": "object",
    }
    print(f"Loading: {sample_path}")
    df = pd.read_csv(sample_path, delimiter='\t', dtype=dtypes)
    df['min_mean_diff'] = df[['mean_diff_father', 'mean_diff_mother']].min(axis=1)

    return df

def load_sample(sample_path: str, min_denovo_coverage: int = 0) -> pd.DataFrame:
    """
    Read and filter TRGT-denovo CSV
    Code source: https://github.com/PacificBiosciences/trgt-denovo/blob/main/scripts/python/trio_analysis.ipynb
    """
    df_sub = read_csv_raw_df(sample_path)
    df_sub = df_sub[df_sub['denovo_coverage'] >= min_denovo_coverage]

    return df_sub


def filter_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter possible de novo expansions
    Code source: https://github.com/PacificBiosciences/trgt-denovo/blob/main/scripts/python/trio_analysis.ipynb
    """
    df_filt = df[
        (df['denovo_coverage'] >= 5) &
        (df['allele_ratio'] >= 0.7) &
        # (df['child_ratio'].between(0.3, 0.7)) &
        (df['mother_dropout_prob'] < 0.005) &
        (df['father_dropout_prob'] < 0.005) &
        (df['mother_dropout'] == "N") & 
        (df['father_dropout'] == "N") &
        (df['mean_diff_mother'] > 7) & 
        (df['mean_diff_father'] > 7)          
    ].copy()        
    df_filt = df_filt[
        (
            (df_filt['child_ratio'].between(0.3, 0.7)) |
            (df_filt['trid'].str.contains('chrX', case=False))
        )
    ]
    df_filt = df_filt.sort_values(by=["min_mean_diff"], ascending=False)

    return df_filt


def compare_to_controls(allele_lens: str, genotype: int, mean: float, std: float, cutoff: float) -> tuple[float, bool]:
    """
    Calculate z score and outlier status for child de novo allele
    """
    def calculate_z_score(value, mean, std):
        try:
            z_score = (value - mean) / std
        except ZeroDivisionError:
            z_score = (value - mean) / 0.0001  # handle division by zero
        return round(z_score, 3)

    denovo_allele_len = get_denovo_al(allele_lens, genotype)

    try:
        z_score = calculate_z_score(denovo_allele_len, float(mean), float(std))
    except: 
        z_score = "."

    outlier = (denovo_allele_len >= float(cutoff)) & (z_score >= 3) # is allele an outlier relative to controls?

    return (round(z_score, 3), outlier)


def get_denovo_al(allele_lens: str, genotype: int) -> int:
    """
    Get length of de novo allele based on genotype, where 1 indicates a de novo short allele and 2 indicates a de novo long allele
    """
    allele_lens = [int(al) for al in allele_lens.split(',')]
    if genotype == 1: # short allele
        denovo_allele_len = min(allele_lens)
    else: # long allele
        denovo_allele_len = max(allele_lens)
    
    return denovo_allele_len


def annotate_segdup(loci: pr.PyRanges, segdup: pr.PyRanges) -> pd.DataFrame:
    """
    Annotate loci against segmental duplications
    """
    loci_segdup = loci.join(
        segdup, suffix="_segdup", how="left", apply_strand_suffix=False
    ).drop(["Strand", "Score"])
    loci_segdup_df = loci_segdup.df
    loci_segdup_df["Segdup"] = np.where(
        loci_segdup_df["Segdup"].isnull(), ".", loci_segdup_df["Segdup"]
    )

    return loci_segdup_df


def calculate_reciprocal_overlap(row):
    """Calculate reciprocal overlap between TRID and constraint region"""
    if pd.isna(row['Start_constraint']) or pd.isna(row['End_constraint']):
        return 0
        
    # Get coordinates
    trid_start = row['Start']
    trid_end = row['End']
    constraint_start = row['Start_constraint'] 
    constraint_end = row['End_constraint']
    
    # Calculate overlap
    overlap_start = max(trid_start, constraint_start)
    overlap_end = min(trid_end, constraint_end)
    overlap_length = overlap_end - overlap_start
    
    # Calculate lengths
    trid_length = trid_end - trid_start
    constraint_length = constraint_end - constraint_start
    
    # Calculate reciprocal overlap
    trid_overlap = overlap_length / trid_length
    try:
        constraint_overlap = overlap_length / constraint_length
    except ZeroDivisionError:
        constraint_overlap = 0
    
    return round(min(trid_overlap, constraint_overlap), 4)

  
def annotate_reg_regions(loci: pd.DataFrame, greendb: str) -> pd.DataFrame:
    """Annotate loci with GREENDB regulatory regions"""
    greendb = pd.read_csv(greendb, sep="\t", compression="gzip")
    greendb.rename(columns={"#Chromosome": "Chromosome"}, inplace=True)
    greendb["Chromosome"] = greendb["Chromosome"].str.replace("chr", "")
    # convert greendb and loci dataframes to PyRanges objects and join them to find overlaps
    greendb_pr = pr.PyRanges(greendb)
    loci_pr = pr.PyRanges(loci)
    loci_pr = loci_pr.join(greendb_pr, how="left", suffix="_greendb")
    loci = loci_pr.df.drop(
        columns=[
            "Start_greendb",
            "End_greendb",
            "regionID",
            "constrain_pct",
            "PhyloP100_median",
            "closestGene_dist",
            "closestProt_symbol",
            "closestProt_dist",
            "N_methods",
        ]
    )
    loci.rename(
        columns={
            "std_type": "GREENDB_reg_region",
            "db_source": "GREENDB_source",
            "closestGene_symbol": "GREENDB_closest_gene",
            "controlled_genes": "GREENDB_controlled_genes",
        },
        inplace=True,
    )
    loci.replace({np.nan: "."}, inplace=True)

    return loci


def group_by_greendb(hits_gene: pd.DataFrame) -> pd.DataFrame:
    """
    One hit may be associated with multiple GREENDB features and therefore multiple rows
    Aggregate by GREENDB and join features to remove duplicate rows
    """
    hits_greendb_dedup = hits_gene.groupby(["trid"]).agg(
        {
            "GREENDB_reg_region": ";".join,
            "GREENDB_source": ";".join,
            "GREENDB_closest_gene": ";".join,
            "GREENDB_controlled_genes": ";".join,
        }
    )
    # merge with original loci table
    hits_gene = hits_gene.drop(
        columns=[
            "GREENDB_reg_region",
            "GREENDB_source",
            "GREENDB_closest_gene",
            "GREENDB_controlled_genes",
        ]
    )
    hits_greendb_merged = hits_gene.merge(hits_greendb_dedup, on=["trid"], how="left")
    hits_greendb_merged_dedup = hits_greendb_merged.drop_duplicates(keep="first")

    return hits_greendb_merged_dedup
