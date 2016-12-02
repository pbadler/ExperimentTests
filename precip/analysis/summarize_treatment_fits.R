rm(list = ls() )

library(stringr)
library(rstan)

# input ------------------------------------------------------------------------------------# 
setwd('~/Documents/ExperimentTests/precip/')
mfiles <- dir('output/stan_fits', 'treatment_effects.*.RDS', full.names = TRUE)

for(i in 1:length(mfiles)){ 
  
  bname <- basename(mfiles[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
  
  temp_fit <- readRDS(mfiles[i])
  
  df <- readRDS(paste0( 'data/temp_data/', vr, '_data_lists_for_stan.RDS'))[[spp]]
  
  model_pars <- temp_fit@model_pars
  
  if(vr == 'growth') { 
    fit_summary <- summary(temp_fit, c('a', 'b1_mu', 'b1', 'bg', 'bt', 'sigma', 'sig_a', 'sig_b1', 'w'))$summary
  }
  
  if(vr == 'survival') { 
    fit_summary <- summary(temp_fit, c('a', 'b1_mu', 'b1', 'bg', 'bt', 'sig_a', 'sig_b1', 'w'))$summary
  }
  
  if(vr == 'recruitment') { 
    fit_summary <- summary(temp_fit, c('a', 'bg', 'bt', 'sig_a', 'theta', 'u', 'w'))$summary
  }
  
  write.csv(fit_summary, paste0( 'output/treatment_model_parameters_', spp, '_', vr, '.csv'))

}

