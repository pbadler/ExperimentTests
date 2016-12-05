rm(list = ls())
library(rstan)
library(stringr)

dl_files <- dir('data/temp_data', 'modified_.*_data_lists_for_stan.RDS', full.names = T)
m_files  <- dir('analysis', 'model.*_1.stan', recursive = T, full.names = T)
best_models <- read.csv('output/best_lppd_scores.csv')

for( i in 1:nrow(best_models)){ 
  
  temp_row <- best_models[i, ]
  m   <- temp_row$model
  spp <- temp_row$species
  vr  <- temp_row$vital_rate
  lambda <- temp_row$lambda
  sd  <- temp_row$sd
  
  dat <- readRDS( dl_files [ grep ( vr, dl_files ) ] ) [[spp]]

  dat$tau_beta <- sd 
  model <- m_files[ grep(vr, m_files )]  
  
  if(is.null( ncol(dat$C)) ){ 
    fit <- stan(str_replace(model, '.stan', '_special.stan'), data  = dat, thin = 4, cores = 4, iter = 2000, seed = 1)
  }else{
    fit <- stan(model, data = dat, thin = 4, cores = 4, iter = 2000, seed = 1)
  }
  
  saveRDS(fit, paste0('output/stan_fits/regularization/', spp[j], '_', vr, '_best_climate_fit.RDS'))
  
  rm( fit )
}

