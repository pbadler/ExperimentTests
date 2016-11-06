#!/bin/bash
#SBATCH --job-name=initialize_models
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --output=SLURM.out
#SBATCH --mail-user=arklein@aggiemail.usu.edu
#SBATCH --constraint=pandorina

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH "--args /projects/A01633220/precip/ data/temp_data/short_model_table_oos.csv $SLURM_ARRAY_TASK_ID 0 $1" run_stan_models_simple_oos.R
