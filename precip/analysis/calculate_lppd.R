rm(list = ls() )

library(stringr)
library(rstan)

# input ------------------------------------------------------------------------------------# 
setwd('~/Documents/ExperimentTests/precip/')

mfiles <- dir('output/stan_fits', 'fit.RDS', full.names = TRUE)
mfiles <- mfiles [ -grep('treatment', mfiles)] # don't use year effects

# log-pointwise predictive density -------------------------------------------------------# 

compute_lppd <- function( stan_fit, ll = 'log_lik' ) { 
  log_lik <- rstan::extract(stan_fit, ll)[[ll]]
  lppd <- log(colMeans(exp(log_lik)))
  lppd
} 

# ---------------------------------------------------------------------------------------------------------------------
i = 2

for( i in 1:length(mfiles)){ 
  
  bname <- basename(mfiles[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
  m <- mpars[3]
  
  temp_fit <- readRDS(mfiles[i])
  dat      <- readRDS(paste0('data/temp_data/modified_', vr, '_data_lists_for_stan.RDS'))[[spp]]
  
  # log-pointwise predictive density ------------------------------------------------------------------------------------# 
  #lppd1 <- compute_lppd(temp_fit)
  lppd2 <- compute_lppd(temp_fit, 'log_lik2')  
  
  rm(temp_fit)
  
  if( vr == 'recruitment'){ 
    y_out <- data.frame(species = spp, vital_rate = vr, model = m, quad = dat$quadhold, trackid = NA, 
                        size = NA, obs_id = dat$obs_idhold, year = dat$yearhold, 
                        Treatment = dat$treathold, lppd2 = lppd2 , Xcenter = NA, Xscale = NA )
  }else{ 
    y_out <- data.frame(species = spp, vital_rate = vr, model = m, quad = dat$quadhold, trackid = dat$trackidhold, 
                        size = dat$Xhold, obs_id = dat$obs_idhold , year = dat$yearhold, 
                        Treatment = dat$treathold, lppd2 = lppd2)
  }
  
  y_out$Treatment <- factor(y_out$Treatment, labels = c('Control', 'Drought', 'Irrigation'))
  
  if( i == 1 ) { 
    write.table( y_out, file = file.path('output', 'lppd_scores.csv'), sep = ',', row.names = FALSE, append = FALSE )
  }else { 
    write.table( y_out, file = file.path('output', 'lppd_scores.csv'), sep = ',', col.names = FALSE, row.names = FALSE, append = TRUE )
  }
  rm(y_out, lppd2, dat)
}

