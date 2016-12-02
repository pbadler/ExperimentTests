rm(list = ls())

library(rstan)

find_dv_trans <- function(x){ 
  ss <-  get_sampler_params(readRDS(x)) 
  
  dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
  
  return(data.frame( model = basename(x), divergent_trans = dv ))
}


myfits <- dir('output/stan_fits', 'fit.RDS', full.names = TRUE)

out <- lapply( myfits, find_dv_trans)

check_df <- do.call(rbind, out)

write.csv(check_df , 'output/check_divergent_transitions.csv')
