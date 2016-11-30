rm(list = ls())
library(rstan)

df <- expand.grid(species = c('ARTR', 'HECO', 'POSE', 'PSSP'), vital_rate = c('growth', 'recruitment', 'survival'))

i = 1

for(i in 1:nrow(df)){ 
  
  spp <- df$species[i]
  vr  <- df$vital_rate[i]
  
  test <- readRDS(paste0('data/temp_data/', vr, '_data_lists_for_stan.RDS'))[[spp]]
  
  myfit <- stan(paste0('analysis/', vr, '/model_', vr, '_treatment_effects.stan'), data = test, cores = 4, iter = 2000, thin = 2,
                control = list(adapt_delta = 0.95, stepsize = 0.5, max_treedepth = 20))
  
  saveRDS(myfit, paste0('output/stan_fits/treatment_effects_', spp, '_', vr, '_.RDS'))
}


