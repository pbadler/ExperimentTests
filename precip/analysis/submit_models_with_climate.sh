#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1-00:00:00
#SBATCH --constraint=gonium
#SBATCH --mail-user=arklein@aggiemail.usu.edu

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH "--args /pscratch/A01633220/precip/ $1 1 2 $SLURM_ARRAY_TASK_ID 1 200 log_lik" run_stan_models.R
