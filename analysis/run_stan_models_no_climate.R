#!/usr/bin/env Rscript

##########################################################################
#
##  Run from the command line given argument for species, vital rate and 
##  number of chains. 
##
##  Runs separate models for each species and vital rate. Only runs the models that 
##  do not involve climate effects, i.e. do not require regularization. 
##
##  Saved output:
##    Save full model output from each model  
##
# 
##########################################################################

rm(list = ls())
library(rstan)

args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==6){
  
  # ---Set working directory, species, vital rate, model number, and number of chains -----------------------------#
  args <- commandArgs(trailingOnly = TRUE)
  
  # example:  args <- c("~/Documents/precip_experiment/", 1, 1, 1, 2)
  
  setwd(args[1])  # set to directory with the "data", "analysis" and "output" folders '/projects/A01633220/precip_experiment/'
  
  do_spp <- args[2] 
  do_vr  <- args[3] 
  do_model <- as.numeric(args[4] )
  nchains <- as.numeric(args[5] )
  niter <- as.numeric(args[6] )

}else if(length(args) == 0 ){ 
  
  # defaults ------------------------------------------------------------------# 
  do_spp <- c('ARTR', 'HECO', 'POSE', 'PSSP')
  do_vr <- c('growth')
  do_model <- c(1, 3)
  nchains <- 4 
  niter <- 2000
  # ---------------------------------------------------------------------------#
  
}else { stop('Incorrect number of arguments supplied.  Supply either 0 or 6.')}


source('analysis/run_stan_model_fxns.R')

for ( i in do_spp ) { 
  for ( m in do_model ) {
    run_stan_model(i, 'growth', m, nchains = nchains, niter = niter)    
  }
}



