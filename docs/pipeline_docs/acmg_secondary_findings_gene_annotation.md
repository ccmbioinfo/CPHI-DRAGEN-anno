## ACMG Secondary Findings Annotation

Michelle Spivak, Madeline Couse

**Version 2026-04**

Variants in CH reports are annotated against an ACMG Secondary Finding gene list tsv. 
This annotation is applied after compound heterozygous status is determined and runs on all CH reports. To run annotation, set `config["run"]["acmg_sf"]` to `true`. The current list being used is ACMG SF v3.3 downloaded from [ClinGen](https://search.clinicalgenome.org/kb/genes/acmgsf?page=1&size=25&order=asc&sort=symbol&search=) 

### Gene matching

Each variant's gene field is matched against ACMG SF gene symbols. Gene fields may contain multiple symbols and a single variant may match multiple ACMG SF genes (e.g fusions or very large CNVs).

### Inputs
- CH reports: wgs.coding, wgs.high.impact, sv, cnv, and wgs.denovo if available. 
- `config["annotation"]["general"]["acmg_sf_list"]`: absolute path to the ACMG SF gene list TSV.
- `config["annotation"]["general"]["acmg_sf_version"]`: version string used in the output column name.

### Output
An `ACMG_SF_{version}` column is appended to each CH report CSV and a .SF suffix is appended to the file name. Variants in genes on the ACMG SF list are flagged with the matching gene symbol(s) separated by `;`. All other variants receive `.`.
A combined family-level ACMG SF report is also generated: `reports/{family}.acmg_sf_report.csv`

This report contains ACMG-flagged variants pulled from: wgs.coding, sv, and cnv reports. High-impact and denovo status is added as a flag when matching variants are also present in the high-impact or denovo CH reports.


### Updating to a new ACMG SF version
Update the TSV path and version string in `config["annotation"]["general"]["acmg_sf_list"]` / `["acmg_sf_version"]`.

## Column descriptions

ACMG SF annotation column in wgs.coding, wgs.high.impact, sv, cnv, and wgs.denovo if available. 

| Column | Description | Example |
| --- | --- | --- |
| `ACMG_SF_v{version}` | ACMG Secondary Findings gene symbol match. Multiple matched genes are separated by `;`. Variants without a match receive `.`. Checked against the ACMG SF gene list TSV provided in `config["annotation"]["general"]["acmg_sf_list"]` | `BRCA1`; `.` |

`reports/{family}.acmg_sf_report.csv`
Combined ACMG SF report columns

| Column | Description | Example |
| --- | --- | --- |
| `position` | Variant position. For SNVs this is copied from the CH report `Position` column. For SV/CNV records this is generated from `CHROM:POS` | `chr9:7542000` |
|`end`| End coordinate for SV/CNV records. SNV records receive `.`| `75420999`|
| `ref` | Reference allele for SNV records. SV/CNV records receive `.` | `G` |
| `alt`| Alternate allele for SNV records. SV/CNV records receive `. `A` |
| `svtype` | Structural variant type for SV/CNV records. SNV records receive `.` | `DEL` |
| `gene`| Gene name or gene symbols associated with the variant | `BRCA1` |
| `acmg_sf_gene` | ACMG SF gene symbol or symbols matched for this variant. Multiple matches are separated by `;`. | `BRCA1` |
| `consequence` | Variant consequence or variant class, depending on source report | `missense_variant`; `DEL` |
| `family` | Family or project identifier | `FAM-001267` |
| `sample_zygosities` | Collapsed zygosity values across available sample zygosity columns. Values are formatted as `sample=value` and separated by `;`. | `CA1C_000001_EXP_0001=Het;CA1C_000002_EXP_0001=Hom` |
| `sample_genotypes` | Collapsed genotype values across available sample genotype columns. Values are formatted as `sample=value` and separated by `;` |`CS1C_000001_EXP_0001=0/1; |
| `clinvar` | ClinVar annotation when available | `Pathogenic` |
| `gnomad_af` | gnomAD allele frequency when available | `0.00001` |
| `ucsc_link` | Link to the locus in the UCSC genome browser| `=HYPERLINK("http://genome.ucsc.edu/...", "UCSC link")` |
| `in_high_impact_report` | Indicates whether the same ACMG-flagged variant was also found in the WGS high-impact report. | `True`; `.`; |
| `variant_reported_in` | Source report where the variant row was pulled from.  | `wgs.coding`; `sv`; `cnv` |