rm(list = ls())
library(rstan)

df <- expand.grid(species = c('ARTR', 'HECO', 'POSE', 'PSSP'), vital_rate = c('growth', 'recruitment', 'survival'))
i = 1
for(i in 1:nrow(df)){ 
  
  spp <- df$species[i]
  vr  <- df$vital_rate[i]
  
  dat <- readRDS(paste0('data/temp_data/modified_', vr, '_data_lists_for_stan.RDS'))[[spp]]
  
  myfit <- stan(paste0('analysis/', vr, '/model_', vr, '_treatment_effects.stan'), data = dat, cores = 4, iter = 2000, thin = 4, seed = 1)
  
  ss <-  get_sampler_params(myfit) 
  
  dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
  
  if ( dv > 0 ) { 
    # try again if divergence 
    inits <- apply(myfit, 2, relist, skeleton = rstan:::create_skeleton(myfit@model_pars, myfit@par_dims)) # find initial values from end of chain 
    myfit <- stan( fit = myfit, init = inits, data = dat, cores = 4, iter = 4000, thin = 8, seed = 1)
  }
  
  ss <-  get_sampler_params(myfit) 
  
  dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
  
  if ( dv > 0 ){ 
    # try again if divergent 
    inits <- apply(myfit, 2, relist, skeleton = rstan:::create_skeleton(myfit@model_pars, myfit@par_dims))
    myfit <- stan( fit = myfit, init = inits, data = dat, cores = 4, iter = 2000, thin = 4, 
                  control = list(adapt_delta = 0.85, stepsize = 0.8, max_treedepth = 20), seed = 1 )
  } 
  
  ss <-  get_sampler_params(myfit) 
  
  dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
  
  if ( dv > 0 ){ 
    # try again if divergent 
    inits <- apply(myfit, 2, relist, skeleton = rstan:::create_skeleton(myfit@model_pars, myfit@par_dims)) 
    myfit <- stan( fit = myfit, init = inits, data = dat, cores = 4, iter = 4000, thin = 8, 
                   control = list(adapt_delta = 0.99, stepsize = 0.1, max_treedepth = 20), seed = 1 )
  } 
  
  saveRDS(myfit, paste0('output/stan_fits/', spp, '_', vr, '_treatment_fit.RDS'))
  
  rm(myfit)

}

