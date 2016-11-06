#!/bin/bash
#SBATCH --job-name=gather_WAIC
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4-00:00:00
#SBATCH --mail-user=arklein@aggiemail.usu.edu

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH gather_WAIC_scores_oos.R
