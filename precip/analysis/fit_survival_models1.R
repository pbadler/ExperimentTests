rm(list = ls())

library(rstan)
library(tidyverse)
library(loo)

source('analysis/stan_data_functions.R')

vr <- 'survival'
mod <- rstan::stan_model('analysis/survival/logistic.stan') # load survival model 

# STAN pars -------------- 
ncores <- 4 
niter <- 2000
nchains <- 4 
nthin <- 1
# --------------------------

# pars used for all species ----------------
k <- 10                    # number of folds 

small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formX = as.formula(paste0 ('~ size + small + W.intra + W.inter + C')) ### Fixed effects design matrix (include climate as "C")
# ------------------------------------------

# set up climate variable table --------------------------# 
load('data/temp_data/climate_combos.RData')

nvars <- length(T_combos)

climate_effects <- paste0( 'C.T.', 1:nvars, '*', 'C.VWC.', 1:nvars)
climate_effects <- c('NULL', climate_effects)

T_combos <- lapply( T_combos , paste, collapse = ', ')
VWC_combos <- lapply( VWC_combos, paste, collapse = ', ')

model_scores <- data.frame( climate_effects, Temperature = c(NA, unlist( T_combos )) , VWC = c(NA, unlist( VWC_combos ))) 
model_scores$out_of_sample_lppd <- NA
model_scores$out_of_sample_mse <- NA
model_scores$ndiv <- NA

total <- nrow( expand.grid( 1:nrow(model_scores), 1:k ) )
# --------------------------------------------------------- #
species <- c('ARTR', 'HECO', 'POSE', 'PSSP')

adapt_delta <- c(0.98, 0.98, 0.8, 0.8)
formX  <- list(formX, formX, formX, formX )

for( s in 1:length(species)){ 

  # choose species 
  sp <- species[s]
  ad <- adapt_delta[s]
  fx <- formX[[s]]

  # Fit survival model -------------------------------------------- # 
  dat_file <- paste0('data/temp_data/', sp, '_growth_survival_dataframe.RDS')
  dat <- readRDS(dat_file)

  dat <- dat[ dat$Period == 'Historical',  ] 
  
  folds <- kfold_split_stratified(k, dat$yid) 
  
  dat$folds <- folds
  
  k_folds <- 
    dat %>% 
    distinct(yid, folds)
  
  dat$size <- scale( dat$logarea.t0 )
  
  intra_comp <- paste0('W.', sp)
  
  dat$size <- scale( dat$logarea.t0 )
  hist(dat$size) ## histogram of size to choose "small" size cutoff 
  
  dat$small <- as.numeric(dat$size < small)
  dat$Y    <- scale( dat$logarea.t1 )
  dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 
  dat$W.intra  <- scale( dat[ , intra_comp])
  dat$W.inter <- scale( rowSums(dat$W[, -( grep ( intra_comp , colnames(dat$W))) ] ) ) # inter specific comp. 
  
  counter <- 1
  
  temp_beta <- list()
  beta_est <- list()
  
  for( j in 1:nrow(model_scores)) { 
    
    # get climate effects 
    formC <- as.formula( paste0 ( '~-1 + ', climate_effects[j]  ))  ### Climate effects design matrix 
    
    lpd <- list()
    mse <- list()
    div <- 0
    
    for( i in 1:k ){
      hold <- k_folds$yid[ k_folds$folds == i  ] 
    
      dl <- process_data(dat = dat, 
                         formX = fx, 
                         formC = formC, 
                         formZ = formZ, 
                         vr = vr, 
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
                               pars = c('hold_log_lik', 'hold_fixef', 'beta'), 
                               control = list(adapt_delta = ad), 
                               refresh = -1)

      div <- div + find_dv_trans(fit1)
      lpd[[i]] <- get_lpd(fit1)
      fixef <- as.numeric( summary(fit1, 'hold_fixef')$summary[,1] )
      mse[[i]] <- mean((inv_logit(fixef) - dl$hold_S)^2)
      
      temp_beta[[i]] <- summary(fit1, 'beta')$summary[,1]
      
      counter <- counter + 1 
      
    }
    
    beta_est[[j]] <- colMeans( do.call(rbind, temp_beta ) )
    model_scores$ndiv[j] <- div
    model_scores$out_of_sample_lppd[j] <- sum( unlist( lpd ) )
    model_scores$out_of_sample_mse[j] <- mean(unlist(mse))
    
  }
  
  beta_est <- do.call(bind_rows, beta_est)
  
  model_scores$spp <- sp
  model_scores$vr <- vr
  
  model_scores <- cbind( model_scores[, 1:8], beta_est)
  
  saveRDS(model_scores, paste0( 'output/', sp, '_', vr, '_model_scores2.RDS'))
  
  model_scores[,4:ncol(model_scores)] <- NA

}

