###########################################################################################
#
#   Read best model   
#   get data  
#   apply prior stdev
#   run model 
#
###########################################################################################

rm(list = ls())

library(dplyr)
library(tidyr)
library(rstan)
library(parallel)
library(gridExtra)


# functions ------------------------------------------------------------------------------------


# input ------------------------------------------------------------------------------------

best_priors <- readRDS('output/best_priors.RDS')

# generate predictions -----------------------------------------------------------------

for ( i in 1:nrow(best_priors)) {
  
  temp_inits <- readRDS(as.character( best_priors$fn[i]))@inits

  do_species <-  best_priors$species[i]
  do_vr <- best_priors$vital_rate[i]
  do_prior_stdev <- best_priors$best_prior[i]
  do_model <- best_priors$model[i]
  
  data_list <- readRDS(file.path('data', 'temp_data', paste(do_vr, 'data_lists_for_stan.RDS', sep = '_')))[[ do_species ]]
  
  data_list$tau_beta <- do_prior_stdev
  
  # select only intraspecific w's for models with intraspecific crowding only 

  if( 'w' %in% names( temp_inits[[1]] )  ) { 
    if( length(temp_inits[[1]]$w) == 1) { 
      w_names <- colnames(data_list$W)
      data_list$W <- data_list$W[ ,  grep(pattern = do_species, w_names) ] 
      data_list$Whold <- data_list$Whold[ , grep(pattern = do_species, w_names) ]
      data_list$W_covs <- 1 
    } 
  }
  
  # fit -------------------------------------------------------------------------------------------# 
  if (do_model %in% c(4,5)) { 
  
  temp_fit <- 
    stan( file = file.path('analysis', do_vr, paste('model', do_vr, do_model , 'predict.stan', sep = '_')), 
        data =  data_list, 
        init = temp_inits, 
        pars = c('a_mu', 'b1_mu', 'b2', 'w', 'mu', 'muhat', 'sigmahat', 'log_lik', 'y_hat'), 
        iter = 2000, 
        chains = length(temp_inits), cores = length(temp_inits))
  
  }else if(do_model == 3 ){ 
    
    temp_fit <- 
      stan( file = file.path('analysis', do_vr, paste('model', do_vr, do_model , 'predict.stan', sep = '_')), 
            data =  data_list, 
            init = temp_inits, 
            pars = c('a_mu', 'b1_mu', 'w', 'mu', 'muhat', 'sigmahat', 'log_lik', 'y_hat'), 
            iter = 2000, 
            chains = length(temp_inits), cores = length(temp_inits))
    
  }else if(do_model == 2){ 
    
    temp_fit <- 
      stan( file = file.path('analysis', do_vr, paste('model', do_vr, do_model , 'predict.stan', sep = '_')), 
          data =  data_list, 
          init = temp_inits, 
          pars = c('a_mu', 'b1_mu', 'b2', 'mu', 'muhat', 'sigmahat', 'log_lik', 'y_hat'), 
          iter = 2000, 
          chains = length(temp_inits), cores = length(temp_inits))

  }else if(do_model ==1 ){ 
    
    temp_fit <- 
      stan( file = file.path('analysis', do_vr, paste('model', do_vr, do_model , 'predict.stan', sep = '_')), 
          data =  data_list, 
          init = temp_inits, 
          pars = c('a_mu', 'b1_mu', 'mu', 'muhat', 'sigmahat', 'log_lik', 'y_hat'), 
          iter = 2000, 
          chains = length(temp_inits), cores = length(temp_inits))
    
    }
  
  
  # output ---------------------------------------------------------------------------------------#

  saveRDS(temp_fit, file.path('output', 'stan_fits', 'predictions', paste(do_species, do_vr, do_model, 'predictions.RDS', sep = '_') ))
  
}

