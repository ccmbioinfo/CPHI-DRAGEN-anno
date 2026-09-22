#!/bin/bash
#SBATCH --job-name=cphi-dragen-anno
#SBATCH --time=50:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%j.out

source /storage/modules/anaconda/2020.11/etc/profile.d/conda.sh
conda activate /srv/shared/conda_envs/snakemake9.16.2


SF="/srv/shared/pipelines/CPHI-DRAGEN-anno/workflow/Snakefile"
CP="/srv/shared/conda_envs/cphi-dragen-anno-snakemake"
SLURM="/srv/shared/pipelines/CPHI-DRAGEN-anno/slurm-profile/"
CONFIG="config_G4RD.yaml"

export XDG_CACHE_HOME="/srv/shared/pipelines/CPHI-DRAGEN-anno/.cache"
mkdir tmp
export TMPDIR=`pwd`/tmp

snakemake --use-conda -s ${SF} --conda-prefix ${CP}  --configfile ${CONFIG} --profile ${SLURM} --executor slurm --verbose -p --rerun-incomplete