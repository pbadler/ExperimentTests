#!/bin/bash
#SBATCH --array=1-nfiles
#SBATCH --job-name=get_WAIC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --constraint=gonium
#SBATCH --output=SLURM.out
#SBATCH --mail-user=arklein@aggiemail.usu.edu

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH "--args /pscratch/A01633220/precip/output/stan_fits $SLURM_ARRAY_TASK_ID" batch_WAIC.R
