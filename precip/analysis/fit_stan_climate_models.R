rm(list = ls())
library(rstan)
library(stringr)

dl_files <- dir('data/temp_data', 'modified_.*_data_lists_for_stan.RDS', full.names = T)
m_files  <- dir('analysis', 'model.*_1.stan', recursive = T, full.names = T)

for( i in 1:length(dl_files)){ 
  
  for( j in 1:4){ 
    
    dl <- readRDS(dl_files[i]) 
    
    vr  <- str_extract(dl_files[i], c('growth', 'recruitment', 'survival'))
    vr  <- vr[!is.na(vr)]
    spp <-  names(dl)
    
    dl <- dl[[j]]
    
    myfit <- stan(m_files[i], data  = dl, thin = 4, cores = 4, iter = 2000, seed = 1)
    
    ss <-  get_sampler_params(myfit) 
    
    dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
    
    if ( dv > 0 ) { 
      # try again if divergence 
      myfit <- stan( fit = myfit, data = dl, cores = 4, iter = 4000, thin = 8, seed = 1)
    }
    
    ss <-  get_sampler_params(myfit) 
    
    dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
    
    if ( dv > 0 ){ 
      # try again if divergent 
      myfit <- stan( fit = myfit, data = dl, cores = 4, iter = 2000, thin = 4, 
                     control = list(adapt_delta = 0.85, stepsize = 0.8, max_treedepth = 20), seed = 1 )
    } 
    
    saveRDS(myfit, paste0('output/stan_fits/', spp[j], '_', vr, '_climate_fit.RDS'))
    rm(myfit, ss)
  }

}

