## ACMG Secondary Findings Annotation

Michelle Spivak, Madeline Couse

**Version 2026-04**

Variants in CH reports are annotated against an ACMG Secondary Finding gene list tsv. 
This annotation is applied after compound heterozygous status is determined and runs on all CH reports. To run annotation, set `config["run"]["acmg_sf"]` to `true`. The current list being used is ACMG SF v3.3 (PMID: 40568962) downloaded from [ClinGen](https://search.clinicalgenome.org/kb/genes/acmgsf?page=1&size=25&order=asc&sort=symbol&search=) 

### Gene matching

Each variant's gene field is matched against ACMG SF gene symbols. Gene fields may contain multiple symbols and a single variant may match multiple ACMG SF genes (e.g fusions or very large CNVs). For SVs and CNVs, only variants overlapping the coding sequence of an ACMG gene(s) are considered.

### Inputs
- CH reports: wgs.coding, wgs.high.impact, sv, cnv, and wgs.denovo if available. 
- `config["annotation"]["general"]["acmg_sf_list"]`: absolute path to the ACMG SF gene list TSV.
- `config["annotation"]["general"]["acmg_sf_version"]`: version string used in the output column name.

### Output
An `ACMG_SF_{version}` column is appended to each CH report CSV and a .SF suffix is appended to the file name. Variants in genes on the ACMG SF list are flagged with the matching gene symbol(s) separated by `;`. All other variants receive `.`.
A combined family-level ACMG SF report is also generated: `reports/{family}.ACMG.SF.hg38.csv`

This report contains ACMG-flagged variants pulled from: wgs.coding, sv, and cnv reports. High-impact and denovo status is added as a flag when matching variants are also present in the high-impact or denovo CH reports.


### Updating to a new ACMG SF version
Update the TSV path and version string in `config["annotation"]["general"]["acmg_sf_list"]` / `["acmg_sf_version"]`.

## Column descriptions

ACMG SF annotation column in wgs.coding, wgs.high.impact, sv, cnv, and wgs.denovo if available. 

| Column | Description | Example |
| --- | --- | --- |
| `ACMG_SF_v{version}` | ACMG Secondary Findings gene symbol match. Multiple matched genes are separated by `;`. Variants without a match receive `.`. Checked against the ACMG SF gene list TSV provided in `config["annotation"]["general"]["acmg_sf_list"]` | `BRCA1`; `.` |

`reports/{family}.ACMG.SF.hg38.csv`
Combined ACMG SF report columns

| Column | Description | Example |
| --- | --- | --- |
| `POSITION` | Variant position. For SNVs this is copied from the CH report `Position` column. For SV/CNV records this is generated from `CHROM:POS` | `chr9:7542000` |
|`END`| End coordinate for SV/CNV records. SNV records receive `.`| `75420999`|
| `REF` | Reference allele for SNV records. SV/CNV records receive `.` | `G` |
| `ALT`| Alternate allele for SNV records. SV/CNV records receive `. `A` |
| `SVTYPE` | Structural variant type for SV/CNV records. SNV records receive `.` | `DEL` |
| `GENE`| Gene name or gene symbols associated with the variant | `BRCA1` |
| `ACMG_SF_GENE` | ACMG SF gene symbol or symbols matched for this variant. Multiple matches are separated by `;`. | `BRCA1` |
| `CONSEQUENCE` | Variant consequence | `missense_variant` |
| `FAMILY` | Family or project identifier | `FAM-001267` |
| `SAMPLE_ZYGOSITIES` | Collapsed zygosity values across available sample zygosity columns. Values are formatted as `sample=value` and separated by `;`. | `CA1C_000001_EXP_0001=Het;CA1C_000002_EXP_0001=Hom` |
| `SAMPLE_zyg` | Sample zygosity |
| `SAMPLE_gt` | Sample genotype |
| `CLINVAR` | ClinVar annotation when available | `Pathogenic` |
| `GNOMAD_AF` | gnomAD allele frequency when available | `0.00001` |
| `UCSC_LINK` | Link to the locus in the UCSC genome browser| `=HYPERLINK("http://genome.ucsc.edu/...", "UCSC link")` |
| `IN_HIGH_IMPACT_REPORT` | Indicates whether the same ACMG-flagged variant was also found in the WGS high-impact report. | `True`; `.`; |
| `VARIANT_REPORTED_IN` | Source report where the variant row was pulled from.  | `wgs.coding`; `sv`; `cnv` |