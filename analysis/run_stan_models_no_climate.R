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
if (length(args) != 7){ 
  stop('####### Incorrect number of arguments supplied ####### \n
        ####### Six arguments required:  
        #######  path: path for folder containing analysis, data and output
        #######  do_spp: species name "ARTR" "HECO" "PSSP" or "POSE", 
        #######          multiple species can be run as space separated list 
        #######          e.g. "ARTR HECO" 
        #######  do_vr:  vital rate name, "growth" "survival" or "recruitment"
        #######  do_model:  model number: 1 - 4. Multiple can be run using R 
        #######             notation e.g. "1:4" or "1,2". 
        #######  do_prior:  prior stdev for regularization of climate effects. 
        #######             Multiple can be run using R notation e.g. "1:24" 
        #######  nchains:  number of chains to run: 1-4 
        #######  niter:  number of iterations to run per chain')
}else if (length(args) == 7){
  
  # ---Set working directory, species, vital rate, model number, and number of chains -----------------------------#
  args <- commandArgs(trailingOnly = TRUE)
  
  setwd(args[1])  # set to directory with the "data", "analysis" and "output" folders '/projects/A01633220/precip_experiment/'
  
  do_spp <- strsplit( args[2], ' ' )[[1]]
  do_vr <- strsplit( args[3], ' ' )[[1]]
  
  do_model <- as.numeric(eval(parse (text = strsplit( args[4], ' '))))
  do_prior_sd  <- as.numeric(eval(parse (text = strsplit( args[5], ' '))))
  nchains <- as.numeric(eval(parse (text = strsplit( args[6], ' '))))
  niter <- as.numeric(eval(parse (text = strsplit( args[7], ' '))))

}

source('analysis/run_stan_model_fxns.R')

for ( i in do_spp ) { 
  for (j in do_vr) {
    for ( k in do_model ) {
      for( l in do_prior_sd ) { 
        run_stan_model(i, j, k, l, nchains = nchains, niter = niter)    
      }
    }
  }
}



