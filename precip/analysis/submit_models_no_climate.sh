#!/bin/bash
#SBATCH --array=1-4
#SBATCH --job-name=no_climate
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --constraint=gonium
#SBATCH --output=SLURM.out
#SBATCH --mail-user=arklein@aggiemail.usu.edu

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH "--args /projects/A01633220/precip/ $SLURM_ARRAY_TASK_ID $1 c(1,3) 1 4 2000" get_WAIC.R
