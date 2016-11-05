#!/bin/bash
#SBATCH --job-name=get_WAIC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4-00:00:00
#SBATCH --output=SLURM.out
#SBATCH --mail-user=arklein@aggiemail.usu.edu

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH "--args /projects/A01633220/precip/ $1 $SLURM_ARRAY_TASK_ID" get_WAIC_simple.R "./slurm-out/get_WAIC_simple.Rout$1$SLURM_ARRAY_TASK_ID"
