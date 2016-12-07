rm(list = ls())
library(stringr)
library(rstan)
# fit year effects ------------------------------------- 

dl_files <- dir('data/temp_data', '*_data_lists_for_stan.RDS', full.names = T)
dl_files <- dl_files[ - grep('modified', dl_files) ] 
m_files  <- dir('analysis', 'year_effects.stan', recursive = T, full.names = T)
i  = 1
j = 1
for( i in 1:length(dl_files)){ 
  
  dl  <- readRDS(dl_files[i])
  vr  <- str_extract(dl_files[i], c('growth', 'recruitment', 'survival'))
  vr  <- vr[!is.na(vr)]
  spp <-  names(dl)

  for( j in 1:length(dl)){ 
    
    fit <- stan(m_files[i], data  = dl[[j]], thin = 4, cores = 4, seed = 1, iter = 2000)
    
    ss <-  get_sampler_params(fit) 
    
    dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
    
    if ( dv > 0 ){ 
      rm(fit)
      fit <- stan(mfiles[i], data = dl[[j]], cores = 4, iter = 2000, thin = 4, 
                    control = list(adapt_delta = 0.95, stepsize = 0.4, max_treedepth = 20), seed = 2 )
    } 
    
    saveRDS(fit, paste0('output/stan_fits/', spp[j], '_', vr, '_year_effects_fit.RDS'))
    
    if(vr == 'recruitment'){
      write.csv(summary(fit, c('a', 'sig_a'))$summary, paste0('output/', spp[j], '_', vr, '_year_effects_table.csv'))
    }else{    
      write.csv(summary(fit, c('a', 'b1', 'b1_mu', 'sig_a', 'sig_b1'))$summary, paste0('output/', spp[j], '_', vr, '_year_effects_table.csv'))
    }
    rm(fit)
  }
}
