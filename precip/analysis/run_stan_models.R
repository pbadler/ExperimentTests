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
##    Save full model output from each model  
##
# 
##########################################################################

rm(list = ls())
library(rstan)

# hard code in list of species and vital rate names 

spp_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')
vr_list <- c('growth', 'survival', 'recruitment')

args <- commandArgs(trailingOnly=TRUE)
#args <- c('~/Documents/ExperimentTests/precip/' , 1 , 1, 1, 1, 0, 1, "'log_lik'", TRUE) # for checking 

# test if there is at least one argument: if not, return an error
if (length(args) != 9){ 
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
        #######  niter:  number of iterations to run per chain
        #######  pars:  parameters to save
        #######  predict:  TRUE/FALSE run prediction version of stan model')
}else if (length(args) == 9){
  
  # ---Set working directory, species, vital rate, model number, and number of chains -----------------------------#

  print(paste( c('path', 'spp', 'vital rate', 'models', 'prior sd', 'nchains', 'iter', 'pars', 'predict = ') , args))
  
  setwd(args[1])  # set to directory with the "data", "analysis" and "output" folders '/projects/A01633220/precip_experiment/'
  
  do_spp <- as.numeric(eval(parse(text = strsplit(args[2], ' '))))
  do_vr <- as.numeric(eval(parse(text = strsplit(args[3], ' '))))
  
  do_spp <- spp_list[do_spp]
  do_vr <- vr_list[do_vr]
  
  do_model <- as.numeric(eval(parse (text = strsplit( args[4], ' '))))
  do_prior_sd  <- as.numeric(eval(parse (text = strsplit( args[5], ' '))))
  nchains <- as.numeric(eval(parse (text = strsplit( args[6], ' '))))
  niter <- as.numeric(eval(parse (text = strsplit( args[7], ' '))))
  
  pars <-  args[8] 
  
  predict <- eval(parse(text = args[9]))  

}


# -- Functions section ----------------------------------------------------------------------------------- 
#  
#   Functions to automate running stan models 
#
# -------------------------------------------------------------------------------------------------------- 

source('analysis/run_stan_fxns.R')


for ( i in do_spp ) { 
  for (j in do_vr) {
    for ( k in do_model ) {
      for( l in do_prior_sd ) { 
        
        if (predict) { 
          output_path <- file.path(getwd(),  'output/stan_fits/predictions')
          
          save_file <- file.path( output_path, paste(i, j, k, l, nchains, sep = '_'))
          
          temp_fit <- run_stan_model(i, j, k, l, nchains = nchains, niter = niter, pars = pars, predict = TRUE)
          saveRDS(temp_fit, file = paste0( save_file, '_predict.RDS'))
          
        }else{ 
          output_path <- file.path(getwd(),  'output/stan_fits')
          
          save_file <- file.path( output_path, paste(i, j, k, l, nchains, sep = '_'))
          
          temp_fit <- run_stan_model(i, j, k, l, nchains = nchains, niter = niter, pars = pars, predict = FALSE)
          saveRDS(temp_fit, file = paste0( save_file, '.RDS'))
          
        }
        
        
      }
    }
  }
}



