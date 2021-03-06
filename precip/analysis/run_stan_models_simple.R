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


rm(list = ls() ) 

library(rstan)

args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 6){ 
  stop('####### Incorrect number of arguments supplied ####### \n
       ####### Arguments required:
       #######  working directory 
       #######  table file :  csv file with model combinations 
       #######  line number : 1 - total combination of models in 
       #######  chains: 1 - 4 
       #######  niter:  number of iterations to run per chain
       #######  predictions: TRUE/FALSE')
  
}else if (length(args) == 6){
  
  # ---Set working directory, species, vital rate, model number, and number of chains -----------------------------#
  
  setwd(args[1])  # set to directory with the "data", "analysis" and "output" folders '/projects/A01633220/precip_experiment/'
  
  models <- read.csv(args[2])
  
  do_line <- as.numeric(eval(parse(text = args[3])))
  
  nchains <- as.numeric(eval(parse (text = strsplit( args[4], ' '))))
  niter <- as.numeric(eval(parse (text = strsplit( args[5], ' '))))
  
  predict <- args[6]  
}

source( 'analysis/run_stan_fxns.R')
source( 'analysis/waic_fxns.R')

if ( do_line <= nrow(models)) { 
  
  line <- models[do_line, ]
  
  species <- line$species 
  vital_rate <- line$vital_rate
  model <- line$model
  prior <- line$prior
  pars_list <- eval (parse(text = as.character( line$pars )) )
  print(pars_list)
  nlambda <- line$nlambda
  
}else{ stop('line number is greater than number of models')}

if( predict ) { 
  output_path <- file.path(getwd(),  'output/stan_fits/predictions/')
  ending <- 'predict.RDS'
  use_pars <- pars_list 
  
}else{ 
  output_path <- file.path(getwd(), 'output/stan_fits/')
  ending <- '.RDS'
  use_pars <- 'log_lik'
}


temp_fit <- run_stan_model(species, vital_rate, model, prior, nchains = nchains, niter = niter, pars = use_pars, predict = predict, nlambda = nlambda)

save_file <- file.path( output_path, paste(species, vital_rate, model, prior, nchains, ending, sep = '_'))

saveRDS( temp_fit , file = save_file )

  
  
