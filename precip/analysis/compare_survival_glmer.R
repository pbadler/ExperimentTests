rm(list = ls())

library(lme4)
library(tidyverse)

source('analysis/stan_data_functions.R')

make_model_list <- function(models_df, basic_fixed = c('Y ~ size'), random = '(size|year)'){ 
  
  model_list <- lapply( apply( models_df, 1, function(x) names(models_df)[as.logical(x) ]), paste, collapse = ' + ')
  models <- lapply( model_list, function(x) paste(basic_fixed, x, random, sep = ' + '))
  return(models)
}

fit_model_list <- function(models) { 
  
  m <- list()
  for( i in 1:length(models) ) { 
    m[[i]] <- glmer( models[[i]], data = dat[dat$Period == 'Historical', ], family = 'binomial')
  }
  return(m)
}

rank_models <- function(m) { 
  aics <- lapply(m, MuMIn::AICc )
  model_scores <- data.frame( model = unlist(models), AICC = unlist(aics))
  
  return( model_scores[order(model_scores$AIC), ] )
  
}

# pars used for all species ----------------
vr <- 'survival'
hold <- c(27:30)          ### Choose the hold out years 
formC <- as.formula(~-1)  ### Climate effects design matrix 

models_df <- expand.grid(W.ARTR = c(T,F), 
                         W.HECO = c(T,F), 
                         W.POSE = c(T,F), 
                         W.PSSP = c(T,F), 
                         W.inter = c(T,F), 
                         GroupP2 = c(T,F), 
                         recent = c(T,F))

# ---------------------------------------------------------------------# 
# Fit ARTR survival model -------------------------------------------- # 
dat_file <- 'data/temp_data/ARTR_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

W <- dat$W
W.inter <- scale(rowSums(W[,2:6]))

dat <- 
  dat %>% 
  dplyr::select( - starts_with('W'))

dat$W.inter <- W.inter
dat <- cbind( dat, W )

dat$GroupP2 <- as.numeric ( dat$Group == 'P2' )
dat$recent  <- as.numeric ( dat$Treatment2 == 'Modern')

dat$size <- scale( dat$logarea.t0 )

hist(dat$size) ## histogram of size to choose "small" size cutoff 

# pars ------------------- 
small <- -1               ### Designate a "small" size theshhold 
# ---------------------- 

dat$small <- as.numeric(dat$size < small)
dat$Y    <- dat$survives

temp_mods <- 
  models_df %>% 
  filter( !( (W.HECO|W.POSE|W.PSSP)))


models <- make_model_list(temp_mods, basic_fixed = 'Y ~ size', random = '(size|year)')
m1     <- fit_model_list(models)
rank   <- rank_models(m1)
rank

rank <- rank[order(rank$AICC), ]

rank %>% head

top_mod <- glmer( as.character( rank$model[1] ), data = dat[dat$Period == 'Historical', ], family = 'binomial')
summary(top_mod)

dat1 <- dat[ dat$Period == 'Historical', ]
dat1$resid <- resid( top_mod )

ggplot( dat1, aes( x = Treatment2, y = resid )) + geom_point() + geom_boxplot()
ggplot( dat1, aes( x = Group, y = resid )) + geom_point() + geom_boxplot()

## Fit HECO ------------------------------------------- # 

dat_file <- 'data/temp_data/HECO_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

W <- dat$W
W.inter <- scale(rowSums(W[,c(1,3:6)]))

dat <- 
  dat %>% 
  dplyr::select( - starts_with('W'))

dat$W.inter <- W.inter
dat <- cbind( dat, W )

dat$GroupP2 <- as.numeric ( dat$Group == 'P2' )
dat$recent  <- as.numeric ( dat$Treatment2 == 'Modern')

dat$size <- scale( dat$logarea.t0 )

hist(dat$size) ## histogram of size to choose "small" size cutoff 

# pars ------------------- 
small <- -1               ### Designate a "small" size theshhold 
# ---------------------- 

dat$small <- as.numeric(dat$size < small)
dat$Y    <- dat$survives

temp_mods <- 
  models_df %>% 
  filter( !( (W.ARTR|W.POSE|W.PSSP)))

models <- make_model_list(temp_mods, basic_fixed = 'Y ~ size', random = '(size|year)')
m1     <- fit_model_list(models)
rank   <- rank_models(m1)

rank <- rank[order(rank$AICC), ]

rank %>% head

top_mod <- glmer( as.character( rank$model[1] ), data = dat[dat$Period == 'Historical', ], family = 'binomial')
summary(top_mod)

dat1 <- dat[ dat$Period == 'Historical', ]
dat1$resid <- resid( top_mod )

ggplot( dat1, aes( x = Treatment2, y = resid )) + geom_point() + geom_boxplot()
ggplot( dat1, aes( x = Group, y = resid )) + geom_point() + geom_boxplot()

## Fit POSE ------------------------------------------- # 

dat_file <- 'data/temp_data/POSE_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

W <- dat$W
W.inter <- scale(rowSums(W[,c(1:2,4:6)]))

dat <- 
  dat %>% 
  dplyr::select( - starts_with('W'))

dat$W.inter <- W.inter
dat <- cbind( dat, W )

dat$GroupP2 <- as.numeric ( dat$Group == 'P2' )
dat$recent  <- as.numeric ( dat$Treatment2 == 'Modern')

dat$size <- scale( dat$logarea.t0 )

hist(dat$size) ## histogram of size to choose "small" size cutoff 

# pars ------------------- 
small <- -1               ### Designate a "small" size theshhold 
# ---------------------- 

dat$small <- as.numeric(dat$size < small)
dat$Y    <- dat$survives

temp_mods <- 
  models_df %>% 
  filter( !( (W.ARTR|W.HECO|W.PSSP)))

models <- make_model_list(temp_mods, basic_fixed = 'Y ~ size', random = '(size|year)')
m1     <- fit_model_list(models)
rank   <- rank_models(m1)

rank <- rank[order(rank$AICC), ]

rank %>% head

top_mod <- glmer( as.character( rank$model[1] ), data = dat[dat$Period == 'Historical', ], family = 'binomial')
summary(top_mod)

dat1 <- dat[ dat$Period == 'Historical', ]
dat1$resid <- resid( top_mod )

ggplot( dat1, aes( x = Treatment2, y = resid )) + geom_point() + geom_boxplot()
ggplot( dat1, aes( x = Group, y = resid )) + geom_point() + geom_boxplot()


## Fit PSSP ------------------------------------------- # 

dat_file <- 'data/temp_data/PSSP_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

W <- dat$W
W.inter <- scale(rowSums(W[,c(1:3,5:6)]))

dat <- 
  dat %>% 
  dplyr::select( - starts_with('W'))

dat$W.inter <- W.inter
dat <- cbind( dat, W )

dat$GroupP2 <- as.numeric ( dat$Group == 'P2' )
dat$recent  <- as.numeric ( dat$Treatment2 == 'Modern')

dat$size <- scale( dat$logarea.t0 )

hist(dat$size) ## histogram of size to choose "small" size cutoff 

# pars ------------------- 
small <- -1               ### Designate a "small" size theshhold 
# ---------------------- 

dat$small <- as.numeric(dat$size < small)
dat$Y    <- dat$survives

temp_mods <- 
  models_df %>% 
  filter( !( (W.ARTR|W.HECO|W.POSE)))

models <- make_model_list(temp_mods, basic_fixed = 'Y ~ size', random = '(size|year)')
m1     <- fit_model_list(models)
rank   <- rank_models(m1)

rank <- rank[order(rank$AICC), ]

rank %>% head

top_mod <- glmer( as.character( rank$model[1] ), data = dat[dat$Period == 'Historical', ], family = 'binomial')
summary(top_mod)

dat1 <- dat[ dat$Period == 'Historical', ]
dat1$resid <- resid( top_mod )

ggplot( dat1, aes( x = Treatment2, y = resid )) + geom_point() + geom_boxplot()
ggplot( dat1, aes( x = Group, y = resid )) + geom_point() + geom_boxplot()

ggplot( dat1, aes( x = factor(year), y = resid)) + geom_point() + geom_boxplot()

as.character(rank$model[1])
simple_mod <- glm(Y~size + W.PSSP + W.inter, data = dat[dat$Period == 'Historical', ], family = 'binomial')

dat1$resid2 <- resid(simple_mod)

ggplot( dat1, aes( x = year, y = resid2)) + 
  geom_point() + 
  stat_summary(fun.y = 'mean', geom = 'line',color = 'blue')

