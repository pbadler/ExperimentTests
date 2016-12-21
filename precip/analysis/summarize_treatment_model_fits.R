rm(list = ls() )

library(stringr)
library(rstan)

# input ------------------------------------------------------------------------------------# 
mfiles <- dir('output/stan_fits', '.*treatment_fit.RDS', full.names = TRUE)

for(i in 1:length(mfiles)){ 
  
  bname <- basename(mfiles[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
  
  if(vr == 'growth') { 
    fit_summary <- summary(readRDS(mfiles[i]), c('a', 'b1_mu', 'b1', 'bg', 'bt', 'sigma', 'sig_a', 'sig_b1', 'w'))$summary
  }
  
  if(vr == 'survival') { 
    fit_summary <- summary(readRDS(mfiles[i]), c('a', 'b1_mu', 'b1', 'bg', 'bt', 'sig_a', 'sig_b1', 'w'))$summary
  }
  
  if(vr == 'recruitment') { 
    fit_summary <- summary(readRDS(mfiles[i]), c('a', 'bg', 'bt', 'sig_a', 'theta', 'u', 'w'))$summary
  }
  
  fit_summary <- data.frame(fit_summary)
  fit_summary$vital_rate <- vr 
  fit_summary$species <- spp 
  write.csv(fit_summary, paste0( 'output/treatment_model_parameters_', spp, '_', vr, '.csv'))

}

