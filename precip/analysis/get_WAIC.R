#!/usr/bin/env 

##########################################################################
#
##  Run from the command line given argument for species, vital rate and 
##  number of chains. 
##
##  Runs separate models for each species and vital rate. Only runs the models that 
##  do not involve climate effects, i.e. do not require regularization. 
##
##  Saved output:
##    Save WAIC score for specified models   
##
# 
##########################################################################

rm(list = ls())
library(rstan)

# hard code in list of species and vital rate names 

spp_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')
vr_list <- c('growth', 'survival', 'recruitment')

args <- commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error
if (length(args) != 7){ 
  stop('####### Incorrect number of arguments supplied ####### \n
        ####### Arguments required:  
        #######  path: path for folder containing analysis, data and output
        #######  do_spp: numeric species id ("ARTR" "HECO" "PSSP" or "POSE"), 
        #######          Multiple can be run using R notation e.g. "1:4" or 
        #######          "1,2". 
        #######  do_vr:  numeric vital rate id ("growth" "survival" or "recruitment")
                         Multiple can be run using R notation e.g. "1:4" or 
        #######          "1,2". 
        #######  do_model:  model number: 1 - 4. Multiple can be run using R 
        #######             notation e.g. "1:5" or "1,2". 
        #######  do_prior:  prior stdev for regularization of climate effects. 
        #######             Multiple can be run using R notation e.g. "1:30" 
        #######  nchains:  number of chains to run: 1-4 
        #######  niter:  number of iterations to run per chain')
  
}else if (length(args) == 7){
  
  # ---Set working directory, species, vital rate, model number, and number of chains -----------------------------#
  args <- commandArgs(trailingOnly = TRUE)
  
  print(paste( c('path', 'spp', 'vital rate', 'models', 'prior sd', 'nchains', 'iter', 'pars') , args))
  
  setwd(args[1])  # set to directory with the "data", "analysis" and "output" folders '/projects/A01633220/precip_experiment/'
  
  do_spp <- as.numeric(eval(parse(text = strsplit(args[2], ' '))))
  do_vr <- as.numeric(eval(parse(text = strsplit(args[3], ' '))))
  
  do_spp <- spp_list[do_spp]
  do_vr <- vr_list[do_vr]
  
  do_model <- as.numeric(eval(parse (text = strsplit( args[4], ' '))))
  do_prior_sd  <- as.numeric(eval(parse (text = strsplit( args[5], ' '))))
  nchains <- as.numeric(eval(parse (text = strsplit( args[6], ' '))))
  niter <- as.numeric(eval(parse (text = strsplit( args[7], ' '))))

}


# -- Functions section ----------------------------------------------------------------------------------- 
#  
#   Functions to automate running stan models 
#
# -------------------------------------------------------------------------------------------------------- 

source('analysis/run_stan_fxns.R')
source('analysis/waic_fxns.R')

for ( i in do_spp ) { 
  for (j in do_vr) {
    for ( k in do_model ) {
      for( l in do_prior_sd ) { 
        
        output_path <- file.path(getwd(),  'output/stan_fits/WAIC_scores/')
        
        save_file <- file.path( output_path, paste(i, j, k, l, nchains, 'WAIC.csv', sep = '_'))
        temp_fit <- run_stan_model(i, j, k, l, nchains = nchains, niter = niter, pars = 'log_lik', predict = FALSE)
        
        temp_waic <- waic(temp_fit)
        write.csv(temp_waic, file = save_file, row.names = FALSE)
        
      }
    }
  }
}



