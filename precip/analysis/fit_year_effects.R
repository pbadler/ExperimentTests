rm(list = ls())
library(stringr)
library(rstan)
# fit year effects ------------------------------------- 

dl_files <- dir('data/temp_data', '*_data_lists_for_stan.RDS', full.names = T)
m_files  <- dir('analysis', 'year_effects.stan', recursive = T, full.names = T)

for( i in 1:length(dl_files)){ 
  
  dl  <- readRDS(dl_files[i])
  vr  <- str_extract(dl_files[i], c('growth', 'recruitment', 'survival'))
  vr  <- vr[!is.na(vr)]
  spp <-  names(dl)

  for( j in 1:length(dl)){ 
    
    fit <- stan(m_files[i], data  = dl[[j]], thin = 4, cores = 4, seed = 1)

    saveRDS(fit, paste0('output/stan_fits/', spp[j], '_', vr, '_year_effects_fit.RDS'))
    
    if(vr == 'recruitment'){
      write.csv(summary(fit, c('a', 'sig_a'))$summary, paste0('output/', spp[j], '_', vr, '_year_effects_table.csv'))
    }else{    
      write.csv(summary(fit, c('a', 'b1', 'b1_mu', 'sig_a', 'sig_b1'))$summary, paste0('output/', spp[j], '_', vr, '_year_effects_table.csv'))
    }
  }
}
