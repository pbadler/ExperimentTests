rm(list = ls())
library(stringr)
library(rstan)
# fit year effects ------------------------------------- 

dl_files <- dir('data/temp_data', '*_data_lists_for_stan.RDS', full.names = T)
dl_files <- dl_files[ - grep('modified', dl_files) ] 
m_files  <- dir('analysis', 'year_effects.stan', recursive = T, full.names = T)

grd <- expand.grid(species = c('ARTR', 'HECO', 'POSE', 'PSSP'), vital_rate = c('growth', 'recruitment', 'survival'))
i = 1
for( i in 1:nrow(grd)){ 
  
  vr <- grd$vital_rate[i]
  spp <- grd$species[i]
  
  dl  <- readRDS(dl_files[grep( vr, dl_files)])[[spp]]
  m <- m_files[grep(vr, m_files)]
  myfit <- stan(m, data  = dl, thin = 4, cores = 4, seed = 1, iter = 2000)
    
  ss <-  get_sampler_params(myfit) 
    
  dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
    
  if ( dv > 0 ){ 
    rm(myfit)
    myfit <- stan(fit = myfit, data = dl, cores = 4, iter = 2000, thin = 4, 
                  control = list(adapt_delta = 0.95, stepsize = 0.4, max_treedepth = 20), seed = 2 )
  } 
    
  saveRDS(myfit, paste0('output/stan_fits/', spp, '_', vr, '_year_effects_fit.RDS'))
    
  if(vr == 'recruitment'){
    write.csv(summary(myfit, c('a', 'sig_a'))$summary, paste0('output/', spp, '_', vr, '_year_effects_table.csv'))
  }else{    
    write.csv(summary(myfit, c('a', 'b1', 'b1_mu', 'sig_a', 'sig_b1'))$summary, paste0('output/', spp, '_', vr, '_year_effects_table.csv'))
  }
  rm(myfit)
}


