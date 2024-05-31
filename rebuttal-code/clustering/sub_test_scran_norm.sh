#!/bin/bash
#SBATCH --time=00:00:00
#SBATCH --account=cfang
#SBATCH --mem=500GB
source ~/env/bin/activate
module load singularity

singularity exec ../pipecomp-singularity.sif Rscript code/clustering/test_scran_norm_on_large_datasets.R code/config_files/alternatives.yaml data/large_experiments/E-MTAB-8060/E-MTAB-8060-SCE.RDS
