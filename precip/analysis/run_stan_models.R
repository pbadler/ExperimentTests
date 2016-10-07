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

args <- commandArgs(trailingOnly=TRUE)

args <- c('ARTR', 'survival', '1', '20', '1', 20, pars = "c('a_mu', 'log_lik')", 
             TRUE, 40)

# test if there is at least one argument: if not, return an error
if (length(args) != 9){ 
  stop('####### Incorrect number of arguments supplied ####### \n
        ####### Arguments required:  
        #######  do_spp: numeric species id ("ARTR" "HECO" "PSSP" or "POSE"), 
        #######          Multiple can be run using R notation e.g. "1:4" or 
        #######          "1,2". 
        #######  do_vr:  numeric vital rate id ("growth" "survival" or "recruitment")
                         Multiple can be run using R notation e.g. "1:4" or 
        #######          "1,2". 
        #######  do_model:  model number: 1 - 4. Multiple can be run using R 
        #######             notation e.g. "1:5" or "1,2". 
        #######  do_prior:  prior stdev for regularization of climate effects. 
        #######             Multiple can be run using R notation 
        #######  nchains:  number of chains to run: 1-4 
        #######  niter:  number of iterations to run per chain
        #######  pars:  parameters to save
        #######  predict:  TRUE/FALSE run prediction version of stan model
        #######  regularization:  nlambda')
  
}else if (length(args) == 9){
  
  # ---Set working directory, species, vital rate, model number, and number of chains -----------------------------#

  print(paste( c( 'spp', 'vital rate', 'models', 'prior sd', 'nchains', 'iter', 'pars', 'predict = ' , 'nlambda =') , args))
  
  do_spp <- args[1]
  do_vr <- args[2]

  do_model <- as.numeric( args[3] ) 
  do_prior_sd  <- as.numeric( args[4] ) 
  nchains <- as.numeric( args[5] ) 
  niter <- as.numeric(args[6])
  
  pars <-  parse( text = args[7])
  
  predict <- args[8]
  
  nlambda <- as.numeric( args[9])

}

# -- Functions section ----------------------------------------------------------------------------------- 
#  
#   Functions to automate running stan models 
#
# -------------------------------------------------------------------------------------------------------- 

source('analysis/run_stan_fxns.R')

dl <- readRDS('data/temp_data/recruitment_data_lists_for_stan.RDS')[[do_spp]]

lambda.set <- exp(seq(-5, 15, length=nlambda))
sd_vec <- sqrt(1/lambda.set) # use sd for stan normal distribution 

dl$prior <- sd_vec[do_prior_sd]

fit <- run_stan_model(do_vr = do_vr, do_spp = do_spp, do_model = do_model, do_prior_sd = do_prior_sd, nchains = nchains, niter = niter, pars = pars,predict = predict, nlambda = nlambda)


