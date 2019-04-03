rm(list = ls())
library(rstan)
library(tidyverse)
library(loo)

source('analysis/stan_data_functions.R')

vr <- 'recruitment'

testing <- F
if( testing ){ 
  
  # STAN pars -------------- 
  ncores <- 1 
  niter <- 1000
  nchains <- 1 
  nthin <- 5
  
}else{
  
  # STAN pars -------------- 
  ncores <- 4 
  niter <- 2000
  nchains <- 4 
  nthin <- 5
  
}

ad <- 0.98
spp <- c('ARTR', 'HECO', 'POSE', 'PSSP')

dl <- fit <- list() 
i <- 1
for(i in 1:length(spp)){ 
  sp <- spp[i]
  window <- 'none'

  dat_file <- paste0('data/temp_data/', sp, '_recruitment_dataframe.RDS')
  dat <- readRDS(dat_file)
  
  dat <- 
    dat %>% 
    filter( Period == 'Historical')
  
  dat <- 
    dat %>% 
    mutate( all_cover = cov.ARTR + cov.HECO + cov.POSE + cov.PSSP)
  
  dat$P1 <- dat[ , paste0('cov.', sp) ]
  dat$P2 <- dat[ , paste0('Gcov.', sp)]
  
  dat$P1_inter <- dat$all_cover - dat$P1
  
  dat$P1_inter <- scale( sqrt( dat$P1_inter) )
  dat$P2 <- scale( sqrt( dat$P2) )
  
  if ( window != 'none' ){   
    moist <- paste0( 'C.VWC.', window)
    therm <- paste0( 'C.T.', window )
      
      # get climate effects 
    formC <- as.formula( paste0 ( '~-1 + ', moist, '*', therm  ))  ### Climate effects design matrix 
  }else if( window == 'none'){ 
      
    formC <- as.formula( '~-1')  
  }

  hold <- 0

  formX <- as.formula( ~ Group + P1_inter + P2 + C)
  
  dl[[i]] <- process_recruitment_data(dat = dat, 
                                   formX = formX, 
                                   formC = formC, 
                                   center = T,
                                   hold = hold, 
                                   IBM = 1)
  
  mod <- rstan::stan_model('analysis/recruitment/recruitment_new.stan')
  
  fit[[i]] <- rstan::sampling(mod, 
                       data = dl[[i]], 
                       chains = nchains, 
                       iter = niter, 
                       cores = ncores,
                       thin = nthin,
                       control = list(adapt_delta = ad))
  
}

traceplot(fit[[1]], 'beta')
traceplot(fit[[2]], 'beta')
traceplot(fit[[3]], 'beta')
traceplot(fit[[4]], 'beta')

summary(fit[[1]], 'beta')$summary
summary(fit[[2]], 'beta')$summary
summary(fit[[3]], 'beta')$summary
summary(fit[[4]], 'beta')$summary


for( i in 1:4){ 
  Y_hat <- summary(fit[[i]], 'Y_hat')$summary
  Y <- dl[[i]]$Y
  N <- dl[[i]]$N

  # posterior predictive Check 
  plot(1:N,  Y)
  segments(x0 = 1:N, x1 = 1:N, y0 =  Y_hat[,4], y1 = Y_hat[,8], col = 'blue') 
  
  plot(Y, Y_hat[,6])
  abline(0,1)
  
  print( (sum( Y > Y_hat[, 8]) + sum(Y < Y_hat[, 4]))/N )
}
