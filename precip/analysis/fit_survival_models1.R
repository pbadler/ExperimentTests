rm(list = ls())

library(rstan)
library(tidyverse)

source('analysis/stan_data_functions.R')

mod <- rstan::stan_model('analysis/survival/logistic.stan') # load survival model 

# STAN pars -------------- 
ncores <- 4 
niter <- 2000 
nchains <- 4 
nthin <- 4
# --------------------------

# pars used for all species ----------------
vr <- 'survival'
hold <- c(27:30)          ### Choose the hold out years 
formC <- as.formula(~-1)  ### Climate effects design matrix 
# ------------------------------------------

# Fit ARTR survival model -------------------------------------------- # 

dat_file <- 'data/temp_data/ARTR_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

dat$size <- scale( dat$logarea.t0 )
hist(dat$size) ## histogram of size to choose "small" size cutoff 

# pars ------------------- 
small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formX = as.formula(~ size + small + W + GroupP2 + Treatment2)  ### Fixed effects design matrix (include climate as "C")
# ---------------------- 

dat$small <- as.numeric(dat$size < small)
dat$Y    <- scale( dat$logarea.t1 )
dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 

dl <- process_data(dat = dat, 
                   formX = formX, 
                   formC = formC, 
                   formZ = formZ, 
                   center = T, 
                   vr = 'survival', 
                   hold = hold )

# save the datalist 
saveRDS(dl, '~/Desktop/ARTR_survival_data.RDS')

sfit1 <- rstan::sampling(mod, 
                         data = dl, 
                         chains = nchains, 
                         iter = niter, 
                         cores = ncores, 
                         thin = nthin )

saveRDS(sfit1, '~/Desktop/ARTR_survival_fit1.RDS')

## Fit HECO ------------------------------------------- # 

dat_file <- 'data/temp_data/HECO_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

dat$size <- scale( dat$logarea.t0 )
hist(dat$size) ## histogram of size to choose "small" size cutoff 

# pars ------------------- 
small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formX = as.formula(~ size + small + W + Treatment2)  ### Fixed effects design matrix (include climate as "C")
# ---------------------- 

dat$small <- as.numeric(dat$size < small)
dat$Y    <- scale( dat$logarea.t1 )
dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 

dl <- process_data(dat = dat, 
                   formX = formX, 
                   formC = formC, 
                   formZ = formZ, 
                   center = T, 
                   vr = 'survival', 
                   hold = hold )


# save the datalist 
saveRDS(dl, '~/Desktop/HECO_survival_data.RDS')

sfit1 <- rstan::sampling(mod, 
                         data = dl, 
                         chains = nchains, 
                         iter = niter, 
                         cores = ncores, 
                         thin = nthin )

saveRDS(sfit1, '~/Desktop/HECO_survival_fit1.RDS')

## Fit POSE ------------------------------------------- # 

dat_file <- 'data/temp_data/POSE_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

dat$size <- scale( dat$logarea.t0 )
hist(dat$size) ## histogram of size to choose "small" size cutoff 

# pars ------------------- 
small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formX = as.formula(~ size + W + GroupP2 + Treatment2)  ### Fixed effects design matrix (include climate as "C")
# ---------------------- 

dat$small <- as.numeric(dat$size < small)
dat$Y    <- scale( dat$logarea.t1 )
dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 

dl <- process_data(dat = dat, 
                   formX = formX, 
                   formC = formC, 
                   formZ = formZ, 
                   center = T, 
                   vr = 'survival', 
                   hold = hold)

# save the datalist 
saveRDS(dl, '~/Desktop/POSE_survival_data.RDS')

sfit1 <- rstan::sampling(mod, 
                         data = dl, 
                         chains = nchains, 
                         iter = niter, 
                         cores = ncores, 
                         thin = nthin)

saveRDS(sfit1, '~/Desktop/POSE_survival_fit1.RDS')

## Fit PSSP ------------------------------------------- # 

dat_file <- 'data/temp_data/PSSP_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

dat$size <- scale( dat$logarea.t0 )
hist(dat$size) ## histogram of size to choose "small" size cutoff 

# pars ------------------- 
small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formX = as.formula(~ size + small + W + GroupP2 + Treatment2)  ### Fixed effects design matrix (include climate as "C")
# ---------------------- 

dat$small <- as.numeric(dat$size < small)
dat$Y    <- scale( dat$logarea.t1 )
dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 

dl <- process_data(dat = dat, 
                   formX = formX, 
                   formC = formC, 
                   formZ = formZ, 
                   center = T, 
                   vr = 'survival', 
                   hold = hold)

# save the datalist 
saveRDS(dl, '~/Desktop/PSSP_survival_data.RDS')

sfit1 <- rstan::sampling(mod, 
                         data = dl, 
                         chains = nchains, 
                         iter = niter, 
                         cores = ncores, 
                         thin = nthin )

saveRDS(sfit1, '~/Desktop/PSSP_survival_fit1.RDS')
