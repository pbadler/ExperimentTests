
rm(list = ls())

library(rstan)
library(tidyverse)
library(loo)

source('analysis/stan_data_functions.R')

fls <- dir('data/temp_data/', '*_growth_survival_dataframe.RDS', full.names = T)

df <- lapply( fls, readRDS ) 
df <- lapply( df, function(x) x[x$Period == 'Historical', ])

Y <- lapply( df, function(x) x$logarea.t1 )

make_size_table <- function( Y ) { 
  data.frame( table(Y) ) %>% arrange( Y ) 
}

size_tabs <- lapply( Y, make_size_table)

lapply( size_tabs, function(x) x[ which.max(x$Freq),  ] ) 
exp(-0.6)

size_tabs[[1]] %>% head
size_tabs[[2]] %>% head( 100 )
size_tabs[[3]] %>% head( 100 )
size_tabs[[4]] %>% head( 200 )


# pars used for all species ----------------
vr <- 'growth'
k <- 10 # number of folds 

small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formX = as.formula(paste0 ('~ size + small + W.intra + W.inter + C')) ### Fixed effects design matrix (include climate as "C")
formE = as.formula(~ size)  ### For growth model, size dependent variance design matrix  
# ------------------------------------------

# set up climate variable table --------------------------# 
load('data/temp_data/climate_combos.RData')

climate_effects <- paste0( 'C.T.', 1:7, '*', 'C.VWC.', 1:7)
climate_effects <- c('NULL', climate_effects)

T_combos <- lapply( T_combos , paste, collapse = ', ')
VWC_combos <- lapply( VWC_combos, paste, collapse = ', ')

model_scores <- data.frame( climate_effects, Temperature = c(NA, unlist( T_combos )) , VWC = c(NA, unlist( VWC_combos ))) 
model_scores$out_of_sample_lppd <- NA
model_scores$out_of_sample_mse <- NA

total <- nrow( expand.grid( 1:nrow(model_scores), 1:k ) )
# --------------------------------------------------------- #
species <- c('ARTR', 'HECO', 'POSE', 'PSSP')

adapt_delta <- c(0.95, 0.9, 0.8, 0.8)
formX  <- list(formX, formX, formX, as.formula( '~ size + small + W.intra + W.inter + GroupP2 + C' ) )
left_cut <- c(-1, -1.3, -1.3, -1.3)

s <- 1
j <- 2
i <- 5

for( s in 1:length(species)){ 
  
  # choose species 
  sp <- species[s]
  ad <- adapt_delta[s]
  fx <- formX[[s]]
  lc <- left_cut[s]
  
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
  
  dat$small <- as.numeric(dat$size < small)
  dat$Y    <- scale( dat$logarea.t1 )
  dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 
  dat$W.intra  <- scale( dat[ , intra_comp])
  dat$W.inter <- scale( rowSums(dat$W[, -( grep ( intra_comp , colnames(dat$W))) ] ) ) # inter specific comp. 
  
  dat <- left_censor_df(dat, left_cut = lc)
  
  counter <- 1
  
  temp_beta <- list()
  beta_est <- list()
  div <- list()
  
  for( j in 1 ) {  #:nrow(model_scores)) { 
    
    # get climate effects 
    formC <- as.formula( paste0 ( '~-1 + ', climate_effects[j]  ))  ### Climate effects design matrix 
    
    lpd <- list()
    mse <- list()
    
    for( i in 1:k ){
      hold <- k_folds$yid[ k_folds$folds == i  ] 
      
      dl <- process_data(dat = dat, 
                         formX = fx, 
                         formC = formC, 
                         formZ = formZ, 
                         center = T, 
                         vr = vr, 
                         hold = hold )
      
      hist( dl$Y, breaks = 100 )
      abline(v = (dl$U), col = 'red')
      abline(v = min(dl$Y_obs), col = 'blue')
      mtext( paste0(sp, ' size, fold ', i))
      
      cat('\n\n')
      
      print( paste( '### ---- species ', s, ' out of ', length(species), ' -------- # '))
      print( paste( '### ---- working on rep', counter, 'of', total, ': ', 100*counter/total, '% done ----------###' ))
      
      cat('\n\n')
    }
  }
}

