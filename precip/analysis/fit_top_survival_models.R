# refit top growth models 

rm(list = ls())

library(rstan)
library(tidyverse)
library(loo)

source('analysis/stan_data_functions.R')

vr <- 'survival'

mod <- rstan::stan_model('analysis/survival/logistic.stan') # load stan model 

# STAN pars -------------- 
ncores <- 4 
niter <- 2000 
nchains <- 4 
nthin <- 4
# --------------------------

small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formX = as.formula(paste0 ('~ size + small + W.intra + W.inter + C')) ### Fixed effects design matrix (include climate as "C")
# ------------------------------------------

# set up climate variable table --------------------------# 
stan_mods <- read_csv('~/Dropbox/projects/ExperimentTests/precip/output/stan_model_ranks.csv')

top_growth_mods <- 
  stan_mods %>% 
  group_by( spp, vr ) %>% 
  filter(oos_lppd == max(oos_lppd), vr == 'survival') %>% 
  select( vr, spp, climate_window)

model_list <- expand.grid( 
  spp = unique( top_growth_mods$spp) , 
  vr = 'survival', 
  model = c('top_model', 'NULL_MOD'))

model_list <- 
  model_list %>% 
  left_join(top_growth_mods) %>% 
  mutate( climate_window = ifelse(model == 'NULL_MOD', 'NULL_MOD', climate_window)) 

# --------------------------------------------------------- #
model_list$adapt_delta <- c(0.98, 0.98, 0.8, 0.8)
model_list$formX <- list( formX  )

formXNULL <- update(formX,  ~ . - C)

model_list <- 
  model_list %>% 
  mutate(formX = ifelse( climate_window == "NULL_MOD", list( formXNULL), formX  )) %>% 
  distinct( spp, vr, adapt_delta, climate_window, formX, left_cut )

i <- 1

for(i in 1:nrow(model_list)){ 
  
  # choose species 
  sp <- model_list$spp[i]
  ad <- model_list$adapt_delta[i]
  fx <- model_list$formX[[i]]
  window <- model_list$climate_window[i]
  
  dat_file <- paste0('data/temp_data/', sp, '_growth_survival_dataframe.RDS')
  dat <- readRDS(dat_file)
  
  intra_comp <- paste0('W.', sp)
  
  dat$size <- scale( dat$logarea.t0 )
  dat$small <- as.numeric(dat$size < small)
  dat$Y    <- scale( dat$logarea.t1 )
  dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 
  dat$W.intra  <- scale( dat[ , intra_comp])
  dat$W.inter <- scale( rowSums(dat$W[, -( grep ( intra_comp , colnames(dat$W))) ] ) ) # inter specific comp. 
  
  if ( window != 'NULL_MOD' ){   
    moist <- paste0( 'C.VWC.', window)
    therm <- paste0( 'C.T.', window )
    
    # get climate effects 
    formC <- as.formula( paste0 ( '~-1 + ', moist, '*', therm  ))  ### Climate effects design matrix 
  }else if( window == 'NULL_MOD'){ 
    
    formC <- as.formula( '~-1')  
  }
  
  hold <- unique( dat$yid [ dat$Period == 'Modern'] ) 
  
  dl <- process_data(dat = dat, 
                     formX = fx, 
                     formC = formC,
                     formZ = formZ, 
                     vr = vr, 
                     hold = hold )
    
  print( paste( '### ---- species', sp, '; climate window', window, '--------------------##'))
  print( paste( '### ---- working on model', i, 'of', nrow(model_list),' -------------###' ))
  
  fit <- rstan::sampling(mod, 
                  data = dl, 
                  chains = nchains, 
                  iter = niter, 
                  cores = ncores,
                  thin = nthin, 
                  pars = c('hold_log_lik', 'beta', 'mu', 'hold_mu'), 
                  control = list(adapt_delta = ad), 
                  refresh = -1 )
  
  saveRDS(dl, file = paste0( 'output/stan_fits/', sp, '_', vr, '_model_data.RDS'))
  saveRDS(fit, file = paste0( 'output/stan_fits/', sp, '_', vr, '_model_', window, '_top_model.RDS'))
  
}
