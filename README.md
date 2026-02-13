# CPHI-DRAGEN-anno

## Installation
- Create Snakemake environment: `conda create -c conda-forge -c bioconda -n snakemake9.16.2 snakemake=9.16.2`
- Activate conda environment: `conda activate snakemake9.16.2`
- Install slurm plugin: `pip install snakemake-executor-plugin-slurm`

## File naming
- CPHI identifiers follow the format `XXXX-000000_EXP_0000`, i.e. `project-participant_EXP_experimentnumber`. DRAGEN VCF sample IDs will conform to this format.
- For joint-genotyped data from a family, DRAGEN VCF names will include a family identifier in the format `FAM-000000`. 