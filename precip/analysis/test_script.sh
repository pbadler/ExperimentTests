#!/bin/bash
#SBATCH --array=1-4
#SBATCH --job-name=initialize_models
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --constraint=gonium

R CMD BATCH '--args ~/Documents/ExperimentTests/precip 2 1 1 1 1 1' /home/andy/Documents/ExperimentTests/precip/analysis/run_stan_models.R output.txt
