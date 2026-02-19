# CPHI-DRAGEN-anno

## Installation
- The steps below are not required if running on the SickKids HPC
- Create Snakemake environment: `conda create -c conda-forge -c bioconda -n snakemake9.16.2 snakemake=9.16.2`
- Activate conda environment: `conda activate snakemake9.16.2`
- Install slurm plugin: `pip install snakemake-executor-plugin-slurm`

## File naming
- CPHI identifiers follow the format `XXXX-000000_EXP_0000`, i.e. `project-participant_EXP_experimentnumber`. DRAGEN VCF sample IDs will conform to this format.
- For joint-genotyped data from a family, DRAGEN VCF names will include a family identifier in the format `FAM-000000`. 

## Running the pipeline
1. Clone this repo to your home directory: `git clone https://github.com/ccmbioinfo/CPHI-DRAGEN-anno/`
2. Make a folder in a directory with sufficient space. Copy over the template files `CPHI-DRAGEN-anno/config/config.yaml`, `CPHI-DRAGEN-anno/config/samples.tsv`, `CPHI-DRAGEN-anno/config/units.tsv`, `CPHI-DRAGEN-anno/workflow/CPHI_DRAGEN_anno.sh`. 

```
$ mkdir FAM-000000
$ cp ~/CPHI-DRAGEN-anno/config/config.yaml ~/CPHI-DRAGEN-anno/config/samples.tsv ~/CPHI-DRAGEN-anno/config/units.tsv ~/CPHI-DRAGEN-anno/workflow/CPHI_DRAGEN_anno.sh FAM-000000
$ cd FAM-000000
```

3. Set up pipeline run: 
- reconfigure `samples.tsv` and `units.tsv` to reflect sample names and input files. Note that because several of the inputs are joint-genotyped, one row in the units.tsv file corresponds to one family, not one sample. `units.tsv` must be configured with the path to: the joint-genotyped (or singleton, if no family members were sequenced) DRAGEN `sequence_variant_vcf`, the joint-genotyped (or singleton) DRAGEN `SV_vcf`, the joint-genotyped (or singleton) DRAGEN `CNV_vcf`, the path to the directory containing the DRAGEN repeat VCFs for all family members `repeat_VCF_dir`, and the path to the directory containing the DRAGEN VNTR VCFs for all family members `VNTR_VCF_dir`. The `samples.tsv` file should contain one row per family member, with the `sample` column corresponding to the sample ID.

`units.tsv` example:
```
family	sequence_variant_vcf	SV_vcf	CNV_vcf	repeat_VCF_dir	VNTR_VCF_dir
FAM-000000	/path/to/sequence_variant_vcf	/path/to/SV_vcf	/path/to/CNV_vcf	/path/to/repeat_VCF_dir	/path/to/VNTR_VCF_dir
```
`particpants.tsv` example:
```
sample
XXXX-000000_EXP_0000
```
- Add paths to the HPO term file and pedigree file to config.yaml. 
- Do a dry run: add a `-n` flag to the Snakemake command in `CPHI_DRAGEN_anno.sh`. This will print out the rules that will be run, but not actually run them.
- If the dry run is successful, run the pipeline: `sbatch CPHI_DRAGEN_anno.sh`.
- If you just want to generate a single report, you can specify the report name in the Snakemake command in `CPHI_DRAGEN_anno.sh`. 
