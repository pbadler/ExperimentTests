#!/usr/bin/env 

##########################################################################
#
##  Run from the command line given argument for species, vital rate and 
##  number of chains. 
##
##  Runs separate models for each species and vital rate. Only runs the models that 
##  do not involve climate effects, i.e. do not require regularization. 
##
##  Output:
##    Save model 
# 
##########################################################################


rm(list = ls() ) 

library(rstan)


setwd('/home/andy/Documents/ExperimentTests/precip/')  # set to directory with the "data", "analysis" and "output" folders '/projects/A01633220/precip_experiment/'
models <- read.csv( 'output/best_WAIC_scores.csv')
nchains <- 4
predict <- TRUE  

for( i in 1:nrow( models )){ 
  
  do_line <- i 
  
  source( 'analysis/run_stan_fxns.R')
  source( 'analysis/waic_fxns.R')
  
  if ( do_line <= nrow(models)) { 
  
  line <- models[do_line, ]
  species <- line$species 
  vital_rate <- line$vital_rate
  model <- line$model
  lambda <- line$lambda
  sd <- line$sd
  pars_list <- eval (parse(text = as.character( line$pars )) )
  print(pars_list)
  niter <- line$niter 
  nlambda <- line$nlambda
  
  }else{ stop('line number is greater than number of models')}

  if( predict ) { 
    output_path <- file.path(getwd(),  'output/stan_fits/predictions/')
    ending <- 'predict.RDS'
    use_pars <- pars_list 
    
  }else{ 
    output_path <- file.path(getwd(), 'output/stan_fits/')
    ending <- '.RDS'
    use_pars <- c('log_lik', 'log_lik2')
  }
  
  temp_fit <- run_stan_model(species, vital_rate, model, do_lambda = lambda, do_prior_sd = sd, nchains = nchains, niter = niter, predict = predict, pars = use_pars)
  
  save_file <- file.path( output_path, paste(species, vital_rate, model, lambda, nchains, ending, sep = '_'))
  
  saveRDS( temp_fit , file = save_file )
  
}    
  
