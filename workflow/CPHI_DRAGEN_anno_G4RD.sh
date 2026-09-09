#!/bin/bash
#SBATCH --job-name=cphi-dragen-anno
#SBATCH --time=50:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x-%j.out

micromamba activate /srv/shared/conda_envs/snakemake9.16.2

SF="/srv/shared/pipelines/CPHI-DRAGEN-anno/workflow/Snakefile"
CP="/srv/shared/conda_envs/cphi-dragen-anno-snakemake"
SLURM="/srv/shared/pipelines/CPHI-DRAGEN-anno/slurm-profile/"
CONFIG="config_G4RD.yaml"

export XDG_CACHE_HOME="/srv/shared/pipelines/CPHI-DRAGEN-anno/.cache" # otherwise, seem to have issues with NFS latency
export TMPDIR=`pwd`

snakemake --use-conda -s ${SF} --conda-prefix ${CP}  --configfile ${CONFIG} --profile ${SLURM} --executor slurm --verbose -p --rerun-incomplete