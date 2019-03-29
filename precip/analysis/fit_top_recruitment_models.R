rm(list = ls())

library(rstan)
library(tidyverse)
library(loo)

source('analysis/stan_data_functions.R')

vr <- 'recruitment'

mod <- rstan::stan_model('analysis/recruitment/recruitment.stan') # load stan model 

# STAN pars -------------- 
ncores <- 1 
niter <- 500 
nchains <- 1 
nthin <- 5
# --------------------------

small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ 1)  ### Year effects design matrix 
formX = as.formula(paste0 ('~ 1 + C')) ### Fixed effects design matrix (include climate as "C")
# ------------------------------------------

# set up climate variable table --------------------------# 
stan_mods <- read_csv('~/Dropbox/projects/ExperimentTests/precip/output/model_ranks_new.csv')

top_mods <- 
  stan_mods %>% 
  group_by( spp, vr ) %>% 
  filter(oos_lppd == max(oos_lppd), vr == vr) %>% 
  select( vr, spp, climate_window)

model_list <- expand.grid( 
  spp = unique( top_mods$spp) , 
  vr = vr, 
  model = c('top_model', 'none'))

model_list <- 
  model_list %>% 
  left_join(top_mods) %>% 
  mutate( climate_window = ifelse(model == 'none', 'none', climate_window)) 

# --------------------------------------------------------- #
model_list$adapt_delta <- c(0.98, 0.98, 0.8, 0.8)
model_list$formX <- list( formX  )

formXNULL <- update(formX,  ~ . - C)

model_list <- 
  model_list %>% 
  mutate(formX = ifelse( climate_window == "none", list( formXNULL), formX  )) %>% 
  filter(!is.na(formX)) %>% 
  distinct( spp, vr, adapt_delta, climate_window, formX )

i <- 1

for(i in 1:nrow(model_list)){ 
  
  # choose species 
  sp <- model_list$spp[i]
  ad <- model_list$adapt_delta[i]
  fx <- model_list$formX[[i]]
  window <- model_list$climate_window[i]
  
  dat_file <- paste0('data/temp_data/', sp, '_recruitment_dataframe.RDS')
  dat <- readRDS(dat_file)
  
  dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 
  
  P1.intra <- paste0('cov.', sp)
  P2.intra <- paste0('Gcov.', sp)
  
  dat$P1.intra <- dat[ , P1.intra]
  dat$P2.intra <- dat[ , P2.intra]
  
  dat$P1.inter <- rowSums(dat$parents1[, - ( grep ( P1.intra, colnames(dat$parents1))) ])
  dat$P2.inter <- rowSums(dat$parents2[, - ( grep ( P2.intra, colnames(dat$parents2))) ])
  
  dat$parents1 <- as.matrix( cbind( dat$P1.intra, dat$P1.inter))
  dat$parents2 <- as.matrix( cbind( dat$P2.intra, dat$P2.inter))
  
  if ( window != 'none' ){   
    moist <- paste0( 'C.VWC.', window)
    therm <- paste0( 'C.T.', window )
    
    # get climate effects 
    formC <- as.formula( paste0 ( '~-1 + ', moist, '*', therm  ))  ### Climate effects design matrix 
  }else if( window == 'none'){ 
    
    formC <- as.formula( '~-1')  
  }
  
  hold <- unique( dat$yid [ dat$Period == 'Modern'] ) 
  
  dl <- process_recruitment_data(dat = dat, 
                                 formX = formX, 
                                 formC = formC, 
                                 formZ = formZ,
                                 center = T,
                                 hold = hold )

  print( paste( '### ---- species', sp, '; climate window', window, '--------------------##'))
  print( paste( '### ---- working on model', i, 'of', nrow(model_list),' -------------###' ))
  
  fit1 <- rstan::sampling(mod, 
                          data = dl, 
                          chains = nchains, 
                          iter = niter, 
                          cores = ncores,
                          thin = nthin,
                          pars = c('hold_log_lik', 'hold_SSE', 'beta', 'w'), 
                          control = list(adapt_delta = ad), 
                          refresh = -1)
  
  saveRDS(dl, file = paste0( 'output/stan_fits/', sp, '_', vr, '_', window, '_model_data.RDS'))
  saveRDS(fit, file = paste0( 'output/stan_fits/', sp, '_', vr, '_', window, '_top_model.RDS'))
  
}
