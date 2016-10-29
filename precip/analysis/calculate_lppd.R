rm(list = ls() )

library(stringr)
library(rstan)

# input ------------------------------------------------------------------------------------# 
setwd('~/Documents/ExperimentTests/precip/')

mfiles <- dir('output/stan_fits/predictions', '4_predict.RDS', full.names = TRUE)

# log-pointwise predictive density -------------------------------------------------------# 

compute_lppd <- function( stan_fit, ll = 'log_lik' ) { 
  log_lik <- rstan::extract(stan_fit, ll)[[ll]]
  lppd <- log(colMeans(exp(log_lik)))
  lppd
} 

# ---------------------------------------------------------------------------------------------------------------------

for( i in 1:length(mfiles)){ 
  
  bname <- basename(mfiles[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
  m <- mpars[3]
  lambda <- mpars[4]
  
  temp_fit <- readRDS(mfiles[i])

  # log-pointwise predictive density ------------------------------------------------------------------------------------# 

  lppd1 <- compute_lppd(temp_fit)
  lppd2 <- compute_lppd(temp_fit, 'log_lik2')  
  
  rm(temp_fit)
  
  y_out <- data.frame(species = spp, vital_rate = vr, model = m, lambda = lambda , lppd1 = lppd1 , lppd2 = lppd2 )
  
  if( i == 1 ) { 
    write.table( y_out, file = file.path('output', 'lppd_scores.csv'), sep = ',', row.names = FALSE, append = FALSE )
  }else { 
    write.table( y_out, file = file.path('output', 'lppd_scores.csv'), sep = ',', col.names = FALSE, row.names = FALSE, append = TRUE )
  }
  rm(y_out, lppd1, lppd2)
}


