#!/bin/bash
#SBATCH --array=1-4
#SBATCH --job-name=initialize_models
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --constraint=gonium
#SBATCH --output=SLURM.out
#SBATCH --mail-user=arklein@aggiemail.usu.edu

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH "--args /pscratch/A01633220/precip/ $SLURM_ARRAY_TASK_ID 1 c(1:5) 1 0 1 log_lik" run_stan_models.R
