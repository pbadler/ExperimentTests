#!/bin/bash
#SBATCH --job-name=predictive_models
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --constraint=gonium
#SBATCH --output=SLURM.out
#SBATCH --mail-user=arklein@aggiemail.usu.edu

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH "--args $SLURM_ARRAY_TASK_ID" generate_prediction.R
