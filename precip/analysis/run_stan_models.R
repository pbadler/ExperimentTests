##########################################################################
#
##  R script for growth as run on Utah State University's
##  High Performance Computing system.
## 
##  Adapted from Andrew Tredennick's script 
##
##  Script takes command line argument 'species' and 'vital_rate'.  Runs 
##  m separate models for each species and vital rate. Models that include
##  climate effects are run for each level of n prior standard deviation 
##  for regularization of the climate effects parameters. Script will launch 
##  n*c models:  where 'c' is the the number of models with climate effects 
##  and 'n' the number of prior standard deviations. 
##
##  Saved output:
##    Save full model output from each model  
##
##  Original Author:      Andrew Tredennick
##  Modified By:          Andrew Kleinhesselink 
##  Email:                arklein@aggiemail.usu.edu
##  Date created:         12-6-2015
##  Date modified:        09-08-2016
# 
##########################################################################

rm(list = ls())

library(rstan)
library(matrixStats)
library(loo )

tweak_inits <- function(inits){ 
  
  inits + runif(length(inits), -0.001, 0.001)
  
}

# -- read in growth data lists and initial values ---------------------------------------------------------# 

# ---Set SD Prior from Command Line Arguments -------------------------------------------------------------#

args <- commandArgs(trailingOnly = TRUE)

do_spp <- args[1]
do_vr  <- args[2]
do_sd  <- args[3]
nchains <- args[4]

# Set up regularization grid ------------------------------------------------------------------------------#

n.beta <- 24
sd_vec <- seq(0.1,1.5, length.out = n.beta)

# read in data and models ---------------------------------------------------------------------------------#
do_vr <- 'growth'
do_spp <- 'ARTR'
data_path <- 'data/temp_data'
model_path <- 'analysis'
output_path <- 'analysis/stan_fits'

data_file <- dir(data_path, pattern = paste0(do_vr, '_data_lists_for_stan.RDS'), full.names = TRUE )
init_file <- dir(data_path, pattern = paste0(do_vr, '_init_vals.RDS'), full.names = TRUE)

data_list <- readRDS( data_file )[do_spp][[1]]
init_vals <- readRDS(init_file)[do_spp][[1]]

models <- dir(file.path(model_path, do_vr), '[0-9].stan', full.names = TRUE)

# -- setup prior regularization only for the models with climate ------------------------------------------# 

m <- length(init_vals)
C <- grep( 'b2', lapply ( init_vals, names )) # models with climate 
l <- c(1:5)[ - grep('b2', lapply ( init_vals, names ) ) ] # models without climate 

# run models -----------------------------------------------------------------------------------------------# 

run_stan_iteration <- function(n.beta, stan_fit, do_spp, do_vr, model_num, datalist, inits, output_path = 'analysis', ... ) { 
  
  datalist$tau_beta <- n.beta # set prior 
  
  save_file <- file.path(output_path, paste( do_spp, do_vr, model_num, 'sd', n.beta, sep = '_') )  
  
  stan_fit <- stan( fit = stan_fit, model_name = basename(save_file), data = datalist, ... )
  
  saveRDS( stan_fit, paste0( save_file, '.RDS'))
  
  return( stan_fit )
  
} 

i = 1 

for( i in C ) { 
  
  temp_data <- data_list
  temp_inits <- init_vals[[i]]
  
  # -- select only intraspecific w's for models with intraspecific crowding only --------------------------------# 
  
  if( 'w' %in% names( temp_inits )  ) { 
    if( length(temp_inits$w) == 1) { 
      
      w_names <- colnames(temp_data$W)
      temp_data$W <- temp_data$W[ ,  grep(pattern = do_spp, w_names) ] 
      temp_data$Whold <- temp_data$Whold[ , grep(pattern = do_spp, w_names) ]
      temp_data$W_covs <- 1 
    } 
  }

  # ---------------------------------------------------------------------------------------------------------#  

  pars <- c('log_lik' , 'log_lik2')
  
  temp_inits <- rep( list( temp_inits), 4 )

  temp_inits <- lapply(temp_inits, lapply, tweak_inits )
  
  empty_fit <- stan( file = models[i], data = temp_data, chains = 0, pars = pars )
  
  n.beta <- 1
  nchains <- 2
  
  test <- run_stan_iteration(n.beta, stan_fit = empty_fit, do_spp, do_vr, model_num = i, datalist = temp_data, 
                             output_path, pars = pars, init = temp_inits, chains = nchains, cores = 4)
  
}



# output ---------------------------------------------------------------------------------------------------#


models = 4*3*24 
models

nchains = models*4 
nchains

size_per_chain = 10 

total_size = nchains*size_per_chain

total_size

total_size/1000


