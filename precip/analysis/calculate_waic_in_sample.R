rm(list = ls() )

library(stringr)
library(rstan)

# input ------------------------------------------------------------------------------------# 

mfiles <- dir('output/stan_fits', 'fit.RDS', full.names = TRUE)

# log-pointwise predictive density -------------------------------------------------------# 

source('analysis/waic_fxns.R')

# compute_lppd <- function( stan_fit, ll = 'log_lik' ) { 
#   log_lik <- rstan::extract(stan_fit, ll)[[ll]]
#   lppd <- log(colMeans(exp(log_lik)))
#   lppd
# } 

# ---------------------------------------------------------------------------------------------------------------------
i = 1
for( i in 1:length(mfiles)){ 
  
  bname <- basename(mfiles[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
  m <- mpars[3]
  
  if(m != 'treatment' ) { 
    temp_fit <- readRDS(mfiles[i])
    dat      <- readRDS(paste0('data/temp_data/modified_', vr, '_data_lists_for_stan.RDS'))[[spp]]
    
    out <- waic(temp_fit)
    out$vital_rate <- vr 
    out$species <- spp 
    out$model <- m
    
    if( i == 1 ) { 
      write.table( out, file = file.path('output', 'waic_scores.csv'), sep = ',', row.names = FALSE, append = FALSE )
    }else { 
      write.table( out, file = file.path('output', 'waic_scores.csv'), sep = ',', col.names = FALSE, row.names = FALSE, append = TRUE )
    }
  }
  
}

