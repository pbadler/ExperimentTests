###########################################################################################
#
#   Read best model   
#   get data  
#   apply prior stdev
#   run model 
#
###########################################################################################

rm(list = ls())

library(rstan)

#

setwd('/pscratch/A01633220/precip')

best_priors <- readRDS('output/best_WAIC_scores.RDS')

nchains <- 4 

# 
total_rows <- nrow(best_priors)

args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 1){ 
  stop('####### Incorrect number of arguments supplied ####### \n
       #######  Arguments required:  
       #######  spp: species name ("ARTR" "HECO" "PSSP" or "POSE") 
       #######  vr:  vital rate ("growth" "survival" or "recruitment")
       #######  m:  model number: (1, 2, 3, 4, 5)')
}else if (length(args) == 1){
  
  # ---Set working directory, species, vital rate, model number, and number of chains -----------------------------#
  args <- commandArgs(trailingOnly = TRUE)
  
  do_row <- as.numeric(args[1])
  if( do_row > total_rows ){ 
    stop('row number out of range!')
  } else { 
    best_priors <- best_priors[do_row, ] 
  }
  
}

# generate predictions -----------------------------------------------------------------

do_species <-  best_priors$species
do_vr <- best_priors$vital_rate
do_model <- best_priors$model

temp_inits <- readRDS(file.path('data', 'temp_data', paste(do_vr, 'init_vals.RDS', sep = '_') ))

temp_inits <- rep( temp_inits[[do_species]][do_model], nchains ) 

# regularization based on Gerber et al. 2015 ---------------------------------------------------------------------# 
nlambda <- 30
lambda.set <- exp(seq(-5, 15, length=nlambda))
sd_vec <- sqrt(1/lambda.set) # use sd for stan normal distribution 
# ----------------------------------------------------------------------------------------------------------------# 

do_prior_stdev <- sd_vec [ best_priors$prior ] 

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
      chains = length(temp_inits))

}else if(do_model == 3 ){ 
  
  temp_fit <- 
    stan( file = file.path('analysis', do_vr, paste('model', do_vr, do_model , 'predict.stan', sep = '_')), 
        data =  data_list, 
        init = temp_inits, 
        pars = c('a_mu', 'b1_mu', 'w', 'mu', 'muhat', 'sigmahat', 'log_lik', 'y_hat'), 
        iter = 2000, 
        chains = length(temp_inits))
    
}else if(do_model == 2){ 
    
  temp_fit <- 
    stan( file = file.path('analysis', do_vr, paste('model', do_vr, do_model , 'predict.stan', sep = '_')), 
        data =  data_list, 
        init = temp_inits, 
        pars = c('a_mu', 'b1_mu', 'b2', 'mu', 'muhat', 'sigmahat', 'log_lik', 'y_hat'), 
        iter = 2000, 
        chains = length(temp_inits))

}else if(do_model ==1 ){ 
    
  temp_fit <- 
    stan( file = file.path('analysis', do_vr, paste('model', do_vr, do_model , 'predict.stan', sep = '_')), 
        data =  data_list, 
        init = temp_inits, 
        pars = c('a_mu', 'b1_mu', 'mu', 'muhat', 'sigmahat', 'log_lik', 'y_hat'), 
        iter = 2000, 
        chains = length(temp_inits))
}

# output ---------------------------------------------------------------------------------------#

saveRDS(temp_fit, file.path('output', 'stan_fits', 'predictions', paste(do_species, do_vr, do_model, 'predictions.RDS', sep = '_') ))
  

