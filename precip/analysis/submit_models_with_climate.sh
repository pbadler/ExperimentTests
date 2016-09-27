#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4-00:00:00
#SBATCH --constraint=gonium
#SBATCH --mail-user=arklein@aggiemail.usu.edu

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH "--args /pscratch/A01633220/precip/ $1 $2 $3 $SLURM_ARRAY_TASK_ID 4 2000" get_WAIC.R
