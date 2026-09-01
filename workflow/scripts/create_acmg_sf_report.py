import pandas as pd
import logging
import os

def log_message(*message):
    if message:
        for i in message:
            logging.info(i)
            print(i)

def get_value(row, col, default="."):
    if col in row.index and pd.notna(row[col]) and str(row[col]).strip() != "":
        return row[col]
    return default

def clean_sample_name(sample_name):
    return str(sample_name).replace("-", "_")

def discover_samples(df):
    zygosity_by_sample = {}
    genotype_by_sample = {}

    for col in df.columns:
        if col.startswith("Zygosity."):
            sample = clean_sample_name(col.replace("Zygosity.", ""))
            zygosity_by_sample[sample] = col
        elif col.endswith("_zyg"):
            sample = clean_sample_name(col[:-4])
            zygosity_by_sample[sample] = col
        elif col.endswith("_GT"):
            sample = clean_sample_name(col[:-3])
            genotype_by_sample[sample] = col

    all_samples = sorted(set(zygosity_by_sample) | set(genotype_by_sample))
    return [
        (sample, zygosity_by_sample.get(sample), genotype_by_sample.get(sample))
        for sample in all_samples
    ]


def per_sample_zygosity_genotype_columns(row, sample_columns):
    columns = {}
    for sample, zygosity_col, genotype_col in sample_columns:
        columns[f"{sample}_zyg"] = get_value(row, zygosity_col) if zygosity_col else "."
        columns[f"{sample}_GT"] = get_value(row, genotype_col) if genotype_col else "."
    return columns


def get_acmg_sf_gene_annotations(acmg_gene_string, acmg_sf_list):
    if pd.isna(acmg_gene_string) or str(acmg_gene_string).strip() in ("", "."):
        return ".", "."

    genes = []
    for gene_value in str(acmg_gene_string).split(";"):
        gene = gene_value.strip()
        if gene:
            genes.append(gene)

    # Match the reported gene or genes directly to rows in the ACMG SF list.
    matches = acmg_sf_list[acmg_sf_list["Gene"].isin(genes)]
    phenotypes = matches["Phenotype"].dropna().astype(str).str.strip()
    phenotypes = phenotypes[(phenotypes != "") & (phenotypes != ".")]
    phenotypes = phenotypes.drop_duplicates()

    inheritances = matches["Inheritance"].dropna().astype(str).str.strip()
    inheritances = inheritances[(inheritances != "") & (inheritances != ".")]
    inheritances = inheritances.drop_duplicates()

    return (
        ";".join(phenotypes) if len(phenotypes) > 0 else ".",
        ";".join(inheritances) if len(inheritances) > 0 else ".",
    )


def make_variant_key(row, input_report_type):
    if input_report_type in ["wgs.coding.CH", "wgs.high.impact.CH", "wgs.denovo.CH"]:
        return f"{get_value(row, 'Position')}:{get_value(row, 'Ref')}:{get_value(row, 'Alt')}"

    if input_report_type in ["sv.CH", "cnv.CH"]:
        chrom = str(get_value(row, "CHROM"))
        if chrom != "." and not chrom.startswith("chr"):
            chrom = f"chr{chrom}"

        return f"{chrom}:{get_value(row, 'POS')}:{get_value(row, 'END')}:{get_value(row, 'SVTYPE')}"

    return "."

def make_acmg_sf_report_rows(df, family, input_report_type, acmg_col, acmg_sf_list,):
    acmg_matches = df[df[acmg_col] != "."].copy()
    sample_columns = discover_samples(acmg_matches)
    report_rows = []

    for _, row in acmg_matches.iterrows():
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
            # Copy the RefSeq change already present in the coding report.
            refseq_change = get_value(row, "Refseq_change")

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
            # RefSeq changes are not currently reported for SVs and CNVs.
            refseq_change = "."

        else:
            continue

        acmg_sf_gene = get_value(row, acmg_col)
        acmg_sf_gene_phenotype, acmg_sf_gene_inheritance = (
            get_acmg_sf_gene_annotations(acmg_sf_gene, acmg_sf_list)
        )
        omim_phenotype = get_value(row, "omim_phenotype_all")
        if omim_phenotype == ".":
            omim_phenotype = get_value(row, "omim_phenotype")

        # OMIM comes from the annotated input; the ACMG values come from the
        # secondary-findings gene list loaded above.
        report_row = {
            "POSITION": position,
            "END": end,
            "REF": ref,
            "ALT": alt,
            "SVTYPE": svtype,
            "REFSEQ_CHANGE": refseq_change,
            "GENE": gene,
            "ACMG_SF_GENE": acmg_sf_gene,
            "OMIM_PHENOTYPE": omim_phenotype,
            "ACMG_SF_GENE_PHENOTYPE": acmg_sf_gene_phenotype,
            "ACMG_SF_GENE_INHERITANCE": acmg_sf_gene_inheritance,
            "CONSEQUENCE": consequence,
            "FAMILY": family,
            "CLINVAR": clinvar,
            "GNOMAD_AF": gnomad_af,
            "UCSC_LINK": ucsc_link,
            "IN_HIGH_IMPACT_REPORT": False,
            "IN_DENOVO_REPORT": False,
            "VARIANT_REPORTED_IN": input_report_type,
            "VARIANT_KEY": make_variant_key(row, input_report_type),
        }
        report_row.update(per_sample_zygosity_genotype_columns(row, sample_columns))
        report_rows.append(report_row)

    return pd.DataFrame(report_rows)

def update_flags(acmg_sf_report, df, input_report_type, acmg_col):
    if input_report_type not in ["wgs.high.impact.CH", "wgs.denovo.CH"]:
        return acmg_sf_report

    if len(acmg_sf_report) == 0:
        return acmg_sf_report
    
    acmg_matches = df[df[acmg_col] != "."].copy()

    if len(acmg_matches) == 0:
        log_message(f"No ACMG SF variants in {input_report_type}; no flags updated")
        return acmg_sf_report

    flag_variant_keys = set(
        acmg_matches.apply(lambda row: make_variant_key(row, input_report_type), axis=1)
    )

    if input_report_type == "wgs.high.impact.CH":
        acmg_sf_report["IN_HIGH_IMPACT_REPORT"] = (
            acmg_sf_report["IN_HIGH_IMPACT_REPORT"].astype(bool)
            | acmg_sf_report["VARIANT_KEY"].isin(flag_variant_keys)
        )

    if input_report_type == "wgs.denovo.CH":
        acmg_sf_report["IN_DENOVO_REPORT"] = (
            acmg_sf_report["IN_DENOVO_REPORT"].astype(bool)
            | acmg_sf_report["VARIANT_KEY"].isin(flag_variant_keys)
        )

    return acmg_sf_report


def infer_input_report_type(path):
    basename = os.path.basename(str(path))

    for report_type in [
        "wgs.coding.CH",
        "wgs.high.impact.CH",
        "wgs.denovo.CH",
        "sv.CH",
        "cnv.CH",
    ]:
        if f".{report_type}.SF." in basename:
            return report_type

    raise ValueError(f"Could not infer input report type from path: {path}")


def main(family, input_reports, output_csv, acmg_sf_version, acmg_sf_tsv):
    logfile = f"logs/report/acmg_sf/{family}.acmg_sf_report.log"
    logging.basicConfig(
        filename=logfile,
        filemode="w",
        level=logging.DEBUG,
        format="%(asctime)s:%(message)s",
        datefmt="%Y-%m-%d %H:%M",
    )

    acmg_col = f"ACMG_SF_v{acmg_sf_version}"
    acmg_sf_list = pd.read_csv(acmg_sf_tsv, sep="\t", usecols=["Gene", "Phenotype", "Inheritance"],)
    acmg_sf_list["Gene"] = acmg_sf_list["Gene"].astype("string").str.strip()
    report_parts = []
    flag_inputs = []

    for input_report in input_reports:
        input_report_type = infer_input_report_type(input_report)
        df = pd.read_csv(input_report)

        if acmg_col not in df.columns:
            raise ValueError(f"{input_report} does not contain expected column {acmg_col}")

        if input_report_type in ["wgs.coding.CH", "sv.CH", "cnv.CH"]:
            rows = make_acmg_sf_report_rows(
                df=df,
                family=family,
                input_report_type=input_report_type,
                acmg_col=acmg_col,
                acmg_sf_list=acmg_sf_list,
            )
            if len(rows) > 0:
                report_parts.append(rows)

        elif input_report_type in ["wgs.high.impact.CH", "wgs.denovo.CH"]:
            flag_inputs.append((df, input_report_type))

    if len(report_parts) > 0:
        acmg_report = pd.concat(report_parts, ignore_index=True)
    else:
        acmg_report = pd.DataFrame(columns=[
            "POSITION",
            "END",
            "REF",
            "ALT",
            "SVTYPE",
            "GENE",
            "ACMG_SF_GENE",
            "OMIM_PHENOTYPE",
            "ACMG_SF_GENE_PHENOTYPE",
            "ACMG_SF_GENE_INHERITANCE",
            "REFSEQ_CHANGE",
            "CONSEQUENCE",
            "FAMILY",
            "SAMPLE_ZYGOSITY",
            "SAMPLE_GENOTYPE",
            "CLINVAR",
            "GNOMAD_AF",
            "UCSC_LINK",
            "IN_HIGH_IMPACT_REPORT",
            "IN_DENOVO_REPORT",
            "VARIANT_REPORTED_IN",
            "VARIANT_KEY",
        ])

    for df, input_report_type in flag_inputs:
        acmg_report = update_flags(
            acmg_sf_report=acmg_report,
            df=df,
            input_report_type=input_report_type,
            acmg_col=acmg_col,
        )
    if "VARIANT_KEY" in acmg_report.columns:
        acmg_report = acmg_report.drop(columns=["VARIANT_KEY"])
    
    acmg_report.to_csv(output_csv, index=False)
    log_message(f"{output_csv} created with {len(acmg_report)} rows")


if __name__ == "__main__":
    family = snakemake.wildcards.family
    input_reports = list(snakemake.input.reports)
    output_csv = snakemake.output.report
    acmg_sf_version = snakemake.params.acmg_sf_version
    acmg_sf_tsv = snakemake.input.acmg_sf_list

    main(family, input_reports, output_csv, acmg_sf_version, acmg_sf_tsv)