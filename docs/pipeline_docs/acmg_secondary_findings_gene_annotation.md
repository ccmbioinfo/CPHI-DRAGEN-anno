## ACMG Secondary Findings Annotation

Michelle Spivak, Madeline Couse

**Version 2026-04**

Variants in CH reports are annotated against an ACMG Secondary Finding gene list tsv. 
This annotation is applied after compound heterozygous status is determined and runs on all CH reports. To run annotation, set `config["run"]["acmg_sf"]` to `true`.

### Gene matching

Each variant's gene field is matched against ACMG SF gene symbols. Gene fields may contain multiple symbols and a single variant may match multiple ACMG SF genes (e.g fusions or very large CNVs).

### Inputs
- CH report CSVs: wgs.coding, panel, panel-flank, wgs.high.impact, sv, cnv, and wgs.denovo if available. 
- `config["annotation"]["general"]["acmg_sf_list"]`: absolute path to the ACMG SF gene list TSV.
- `config["annotation"]["general"]["acmg_sf_version"]`: version string used in the output column name.

### Output
An `ACMG_SF_{version}` column is appended to each CH report CSV and a .SF suffix is appended to the file name. Variants in genes on the ACMG SF list are flagged with the matching gene symbol(s) separated by `;`. All other variants receive `.`.

### Updating to a new ACMG SF version
Update the TSV path and version string in `config["annotation"]["general"]["acmg_sf_list"]` / `["acmg_sf_version"]`.