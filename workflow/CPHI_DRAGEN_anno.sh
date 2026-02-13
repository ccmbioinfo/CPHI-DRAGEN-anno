#!/bin/bash
#SBATCH --job-name=cphi-dragen-anno
#SBATCH --time=50:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%j.out

module load anaconda3/4.4.0
conda activate /hpf/largeprojects/ccm_dccforge/dccdipg/Common/conda_envs/snakemake9.16.2


SF=~/CPHI-DRAGEN-ANNO/Snakefile
CP="/hpf/largeprojects/ccm_dccforge/dccdipg/Common/cphi-dragen-anno-snakemake"
SLURM=~/CPHI-DRAGEN-ANNO/slurm_profile/
CONFIG="config.yaml"

snakemake --use-conda -s ${SF} --conda-prefix ${CP}  --configfile ${CONFIG} --profile ${SLURM} --executor slurm --verbose -p --rerun-incomplete
