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

tweak_inits <- function(inits){ 
  
  inits + runif(length(inits), -0.001, 0.001)
  
}

project_wd <- '~/Documents/precip_experiment/' # set to directory with the "data", "analysis" and "output" folders
setwd(project_wd)

data_path <- 'data/temp_data'
model_path <- 'analysis'
output_path <- 'output/stan_fits'

# ---Set species and vital rate ---------------------------------------------------------------------------#

args <- commandArgs(trailingOnly = TRUE)

do_spp <- args[1]
do_vr  <- args[2]
nchains <- args[3]

#do_vr <- 'growth'
#do_spp <- 'ARTR'

# read in data and models ---------------------------------------------------------------------------------#

data_file <- dir(data_path, pattern = paste0(do_vr, '_data_lists_for_stan.RDS'), full.names = TRUE )
init_file <- dir(data_path, pattern = paste0(do_vr, '_init_vals.RDS'), full.names = TRUE)

data_list <- readRDS( data_file )[do_spp][[1]]
init_vals <- readRDS(init_file)[do_spp][[1]]

models <- dir(file.path(model_path, do_vr), '[0-9].stan', full.names = TRUE)

# -- select models without climate effects ----------------------------------------------------------------------# 

m <- length(init_vals)
C <- grep( 'b2', lapply ( init_vals, names )) # models with climate 
l <- c(1:5)[ - grep('b2', lapply ( init_vals, names ) ) ] # models without climate 


# Fit model -----------------------------------------------------------------------------------------------------#

for( i in l ) {           
  
  save_file <- file.path( output_path, paste(do_spp, do_vr, i, sep = '_'))
  
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
  
  
  # -------------------------------------------------------------------------------------------------------------#
  
  temp_fit <- stan(file = models[i], model_name = basename(save_file), init = init_vals[[i]], data = temp_data, chains = nchains, cores = max( 1, nchains) )
  
  # -- output ---------------------------------------------------------------------------------------------------#
  
  saveRDS(temp_fit, file = paste0( save_file, '.RDS'))
}
