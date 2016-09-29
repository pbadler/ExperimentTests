#!/bin/bash
#SBATCH --job-name=predictive_models
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4-00:00:00
#SBATCH --output=SLURM.out
#SBATCH --mail-user=arklein@aggiemail.usu.edu

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH "--args /projects/A01633220/precip/ output/best_WAIC_scores.csv $SLURM_ARRAY_TASK_ID 4 2000 TRUE" run_stan_models_simple.R