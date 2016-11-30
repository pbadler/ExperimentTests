rm(list = ls())
library(rstan)
library(stringr)
# fit year effects ------------------------------------- 

dl_files <- dir('data/temp_data', 'modified_.*_data_lists_for_stan.RDS', full.names = T)
m_files  <- dir('analysis', 'model.*_1.stan', recursive = T, full.names = T)

for( i in 1:length(dl_files)){ 
  
  dl  <- readRDS(dl_files[i])
  
  vr  <- str_extract(dl_files[i], c('growth', 'recruitment', 'survival'))
  vr  <- vr[!is.na(vr)]
  spp <-  names(dl)
  
  for( j in 1:length(dl)){ 
    
    if(is.null( ncol(dl[[j]]$C)) ){ 
      fit <- stan(str_replace(m_files[i], '.stan', '_special.stan'), data  = dl[[j]], thin = 2, cores = 4, iter = 2000)
    }else{
      fit <- stan(m_files[i], data  = dl[[j]], thin = 2, cores = 4, iter = 2000)
    }
    saveRDS(fit, paste0('output/stan_fits/', spp[j], '_', vr, '_climate_fit.RDS'))
    
  }
}

