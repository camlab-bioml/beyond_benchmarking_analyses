#!/bin/bash
#SBATCH --time=00:00:00
#SBATCH --account=cfang
source ~/env/bin/activate
module load singularity

snakemake -s ~/bb-rebuttal/code/rebuttal_snakefile --unlock --jobs 5000 --use-singularity --profile ~/slurm

snakemake -s ~/bb-rebuttal/code/rebuttal_snakefile  --rerun-incomplete  --jobs 5000 -k  --use-singularity --profile ~/slurm
