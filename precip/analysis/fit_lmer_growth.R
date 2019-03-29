rm(list = ls())

library(tidyverse)
library(stringr)
library(lme4)
library(loo)

source('analysis/stan_data_functions.R')

vr <- 'growth'

# pars used for all species ----------------
k <- 10           # number of folds 

small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formX = as.formula(paste0 ('~ size + small + W.intra + W.inter + C')) ### Fixed effects design matrix (include climate as "C")
formE = as.formula(~ size)  ### For growth model, size dependent variance design matrix  
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

basic_form <- 'Y ~ size + small + W.intra + W.inter + '
random_effects <- ' + (size|year)'

total <- nrow( expand.grid( 1:nrow(model_scores), 1:k ) )
# --------------------------------------------------------- #
species <- c('ARTR', 'HECO', 'POSE', 'PSSP')

formX  <- list(formX, formX, formX, formX )
left_cut <- c(-1, -1.3, -1.3, -1.3)

m <- list()

for( s in 1:length(species)){ 
  
  m[[s]] <- list()
  
  # choose species 
  sp <- species[s]
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
  
  intra_comp <- paste0('W.', sp)
  
  dat$size <- scale( dat$logarea.t0 )
  
  dat$small <- as.numeric(dat$size < small)
  dat$Y    <- scale( dat$logarea.t1 )
  dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 
  dat$W.intra  <- scale( dat[ , intra_comp])
  dat$W.inter <- scale( rowSums(dat$W[, -( grep ( intra_comp , colnames(dat$W))) ] ) ) # inter specific comp. 
  
  dat <- left_censor_df(dat, left_cut = lc)
  
  temp_beta <- list()
  beta_est <- list()
  
  counter <- 1
  
  for( j in 1:nrow(model_scores)) { 
    
    # get climate effects 
    formC <- climate_effects[j]  ### Climate effects design matrix 
    

    temp_formula <- paste( basic_form , formC, random_effects  )     
    temp_formula <- str_remove( temp_formula, '\\+  NULL')
    m[[s]][[j]] <- lmer( temp_formula, data = dat)

      
  }
}

model_grid <- expand.grid(species = species, model = c(1:8)) %>% arrange( species, model )
test <- do.call( bind_rows, unlist( rapply(m, fixef, how = 'list'), recursive = F) )

scores <- rapply( m, AIC)
test <- cbind( model_grid, AIC = scores, test )

test <- 
  test %>% 
  gather( par, est, starts_with('C')) %>% 
  filter( !(model > 1 & is.na(est)))

test <- 
  test %>% 
  mutate(Temp = str_detect(par, 'T'), 
         Moist = str_detect(par, 'VWC'), 
         Temp_x_Moist = str_detect(par, '\\:')) %>% 
  mutate( Temp = ifelse( Temp_x_Moist, F, Temp), 
          Moist = ifelse( Temp_x_Moist, F, Moist)) %>%
  gather( type, val, Temp:Temp_x_Moist) %>% 
  filter( val ) %>% 
  select(-val, -par) %>% 
  distinct() %>% 
  spread( type , est)

var_grid <- cbind( model_grid , climate_effects)

test <- 
  test %>% 
    mutate( vr = 'growth') %>% 
    mutate( climate_window = model - 1 ) %>% 
    mutate( climate_window = ifelse( climate_window == 0 , 'NULL_MOD', climate_window)) %>% 
    rename( 'spp' = species, 
            'intercept' = `(Intercept)`, 
            'small_plants' = small, 
            'intra_comp' = W.intra, 
            'inter_comp' = W.inter) %>% 
    arrange( spp, AIC) %>% 
    select( vr, spp, climate_window, AIC, intercept:inter_comp, Temp, Moist, Temp_x_Moist)

write_csv(test, path = '~/Desktop/lmer_growth_model_ranks.csv')
