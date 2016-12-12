rm(list = ls())
library(rstan)
library(stringr)

dl_files <- dir('data/temp_data', 'modified_.*_data_lists_for_stan.RDS', full.names = T)
m_files  <- dir('analysis', 'model.*_1.stan', recursive = T, full.names = T)
j= 1
i= 1
for( i in 1:length(dl_files)){ 
  
  dl  <- readRDS(dl_files[i])
  vr  <- str_extract(dl_files[i], c('growth', 'recruitment', 'survival'))
  vr  <- vr[!is.na(vr)]
  spp <-  names(dl)
  
  for( j in 1:length(dl)){ 
    
    if(is.null( ncol(dl[[j]]$C)) ){ 
      myfit <- stan(str_replace(m_files[i], '.stan', '_special.stan'), data  = dl[[j]], thin = 4, cores = 4, iter = 2000, seed = 1)
    }else{
      myfit <- stan(m_files[i], data  = dl[[j]], thin = 4, cores = 4, iter = 2000, seed = 1)
    }
    
    ss <-  get_sampler_params(myfit) 
    
    dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
    
    if ( dv > 0 ) { 
      # try again if divergence 
      myfit <- stan( fit = myfit, data = dl[[j]], cores = 4, iter = 4000, thin = 8, seed = 1)
    }
    
    ss <-  get_sampler_params(myfit) 
    
    dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
    
    if ( dv > 0 ){ 
      # try again if divergent 
      myfit <- stan( fit = myfit, data = dl[[j]], cores = 4, iter = 2000, thin = 4, 
                     control = list(adapt_delta = 0.85, stepsize = 0.8, max_treedepth = 20), seed = 1 )
    } 
    
    saveRDS(myfit, paste0('output/stan_fits/', spp[j], '_', vr, '_climate_fit.RDS'))
    
  }
  rm( myfit )
}

