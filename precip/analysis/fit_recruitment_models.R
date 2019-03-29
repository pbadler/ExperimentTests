rm(list = ls())

library(rstan)
library(tidyverse)
library(loo)

source('analysis/stan_data_functions.R')

vr <- 'recruitment'
species <- c('ARTR', 'HECO', 'POSE', 'PSSP')

k <- 2                      ### number of folds 

n_mods <- 2
species <- species[1:2]

# STAN pars -------------- 
ncores <- 1 
niter <- 500
nchains <- 1 
nthin <- 5

adapt_delta <- c(0.98, 0.98, 0.8, 0.8)  #  species specific 
# --------------------------

# Model Parameters 
formZ = as.formula(~ 1 )  ### Year effects design matrix 
formX = as.formula(paste0 ('~ 1 + C')) ### Fixed effects design matrix (include climate as "C")

# set up model selection table --------------------------# 
load('data/temp_data/climate_combos.RData')

nvars <- length(T_combos)

climate_effects <- paste0( 'C.T.', 1:nvars, '*', 'C.VWC.', 1:nvars)
climate_effects <- c('NULL', climate_effects)

model_combos <- data.frame( climate_effects = climate_effects)

model_combos <- model_combos %>% head( n_mods )

# --------------------------------------------------------- #

mod <- rstan::stan_model('analysis/recruitment/recruitment.stan') # load stan model 

# ---------------------------------------------------------- # 

total <- k*nrow(model_combos)*length(species)  ### Total number of models to fit 

counter <- 1

for( s in 1:length(species)){ 
  
  sp <- species[s]
  ad <- adapt_delta[s]
  
  dat_file <- paste0('data/temp_data/', sp, '_recruitment_dataframe.RDS')
  
  dat <- readRDS(dat_file)
  dat <- dat[ dat$Period == 'Historical',  ] 
  
  folds <- kfold_split_stratified(k, dat$yid) 
  dat$folds <- folds
  
  k_folds <- 
    dat %>% 
    distinct(yid, folds)
  
  dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 
  
  P1.intra <- paste0('cov.', sp)
  P2.intra <- paste0('Gcov.', sp)

  dat$P1.intra <- dat[ , P1.intra]
  dat$P2.intra <- dat[ , P2.intra]
  
  dat$P1.inter <- rowSums(dat$parents1[, - ( grep ( P1.intra, colnames(dat$parents1))) ])
  dat$P2.inter <- rowSums(dat$parents2[, - ( grep ( P2.intra, colnames(dat$parents2))) ])
  
  dat$parents1 <- as.matrix( cbind( dat$P1.intra, dat$P1.inter))
  dat$parents2 <- as.matrix( cbind( dat$P2.intra, dat$P2.inter))
  
  ## Set up output table 
  temp_scores <- model_combos
  temp_scores$spp <- sp 
  temp_scores$vr <- vr
  temp_scores$climate_window <- as.numeric( str_extract(model_combos$climate_effects, '\\d+'))
  temp_scores$ndiv <- NA
  temp_scores$oos_mse <- NA
  temp_scores$oos_lppd <- NA

  for( j in 1:nrow(temp_scores)) { 
    
    # get climate effects 
    formC <- as.formula( paste0 ( '~-1 + ', model_combos$climate_effects[j]  ))  ### Climate effects design matrix 
    
    lpd <- NA
    sse <- NA
    div <- 0
    
    for( i in 1:k ){
      hold <- k_folds$yid[ k_folds$folds == i  ] 
      
      dl <- process_recruitment_data(dat = dat, 
                         formX = formX, 
                         formC = formC, 
                         formZ = formZ,
                         center = T,
                         hold = hold )
      
      cat('\n\n')
      
      print( paste( '### ---- species ', s, ' out of ', length(species), ' -------- # '))
      print( paste( '### ---- working on rep', counter, 'of', total, ': ', 100*(counter - 1)/total, '% done ----------###' ))
      
      cat('\n\n')
      
      fit1 <- rstan::sampling(mod, 
                              data = dl, 
                              chains = nchains, 
                              iter = niter, 
                              cores = ncores,
                              thin = nthin,
                              pars = c('log_lik', 'hold_log_lik', 'hold_SSE'), 
                              control = list(adapt_delta = ad), 
                              refresh = -1)
      
      div <- div + find_dv_trans(fit1)
      
      lpd[i] <- sum(get_lpd(fit1))
      sse[i] <- summary(fit1, 'hold_SSE')$summary[, 'mean']        
      counter <- counter + 1 
      
    }
    
    temp_scores$ndiv[j] <- div
    temp_scores$oos_lppd[j] <- sum( lpd )
    temp_scores$oos_mse[j]  <- sum(sse)/length(dl$N + dl$hold_N) # divide by total number of observations
    
  }
  
  saveRDS(temp_scores, paste0( 'output/', sp, '_', vr, '_model_scores_new.RDS'))
  
}



