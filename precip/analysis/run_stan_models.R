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


# test if there is at least one argument: if not, return an error
if (length(args) != 8){ 
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
        #######             Multiple can be run using R notation e.g. "1:24" 
        #######  nchains:  number of chains to run: 1-4 
        #######  niter:  number of iterations to run per chain
        #######  pars:  parameters to save')
}else if (length(args) == 8){
  
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
  pars <- as.character(args[8])
  

}


# -- Functions section ----------------------------------------------------------------------------------- 
#  
#   Functions to automate running stan models 
#
# -------------------------------------------------------------------------------------------------------- 

tweak_inits <- function(inits){ 
  # adds random deviation to list of inital values 
  inits + runif(length(inits), -0.001, 0.001)
  
}

run_stan_model <- function(do_spp, do_vr, do_model, do_prior_sd, nchains, niter, pars ) { 
  
  require(rstan)
  
  data_path <- file.path(getwd(), 'data/temp_data')
  model_path <- file.path(getwd(), 'analysis')
  output_path <- file.path(getwd(),  'output/stan_fits')
  
  # read in data and models ---------------------------------------------------------------------------------------#
  
  data_file <- dir(data_path, pattern = paste0(do_vr, '_data_lists_for_stan.RDS'), full.names = TRUE )
  init_file <- dir(data_path, pattern = paste0(do_vr, '_init_vals.RDS'), full.names = TRUE)
  
  initial_fit <- dir(output_path, pattern = paste0(paste(do_spp, do_vr, do_model, 1, 0, sep = '_'), '.RDS'), full.names = TRUE) # check pre-fit models
  
  data_list <- readRDS( data_file )[[do_spp]]
  init_vals <- readRDS(init_file)[[do_spp]]
  
  if(length(data_list) == 0) { stop('No data!!!')}
  if(length(init_vals) == 0) { stop('No init vals!!!')}
  
  models <- dir(file.path(model_path, do_vr), '[0-9].stan', full.names = TRUE)
  
  # -- select model and initial vals -----------------------------------------------------------------------------# 
  
  m <- models [do_model]
  
  temp_inits <- init_vals[[do_model]]
  
  # modify data list for model ------------------------------------------------------------------------------------#
  
  n.beta <- 30
  sd_vec <- seq(0.1,1.5,length.out = n.beta)
  
  data_list$tau_beta <- sd_vec[do_prior_sd]  # modify prior for regularization 
  
  #
  # select only intraspecific w's for models with intraspecific crowding only 
  #
  
  if( 'w' %in% names( temp_inits )  ) { 
    if( length(temp_inits$w) == 1) { 
      w_names <- colnames(data_list$W)
      data_list$W <- data_list$W[ ,  grep(pattern = do_spp, w_names) ] 
      data_list$Whold <- data_list$Whold[ , grep(pattern = do_spp, w_names) ]
      data_list$W_covs <- 1 
    } 
  }
  
  # -- inital values ----------------------------------------------------------------------------------------------#
  
  temp_inits <- rep(list(temp_inits), max(1, nchains) )
  temp_inits <- lapply(temp_inits, lapply, tweak_inits)
  
  # -- run stan ---------------------------------------------------------------------------------------------------#  
  
  save_file <- file.path( output_path, paste(do_spp, do_vr, do_model, do_prior_sd, nchains, sep = '_'))
  
  print(paste('data list has', length(data_list), 'variables'))
  print(paste('init vals has', length(temp_inits[[1]]), 'variables'))
  print(paste('pars requested', pars))
  
  if (length( initial_fit) == 1 ) { 
    print(paste('initial fit being used', initial_fit ))
    initial_fit <- readRDS(initial_fit)
    temp_fit <- stan (fit = initial_fit, model_name = basename(save_file), data = data_list, chains = nchains, init = temp_inits, iter = niter , pars = pars, cores = max(1, nchains))
  }else{   
    temp_fit <- stan (file = m, model_name = basename(save_file), data = data_list, chains = nchains, init = temp_inits, iter = niter , pars = pars, cores = max(1, nchains))
  }
  
  # -- output -----------------------------------------------------------------------------------------------------#
  
  saveRDS(temp_fit, file = paste0( save_file, '.RDS'))
  
}

for ( i in do_spp ) { 
  for (j in do_vr) {
    for ( k in do_model ) {
      for( l in do_prior_sd ) { 
        run_stan_model(i, j, k, l, nchains = nchains, niter = niter, pars = pars)    
      }
    }
  }
}



