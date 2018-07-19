rm(list = ls())

library(rstan)
library(tidyverse)

source('analysis/stan_data_functions.R')

gmod <- rstan::stan_model('analysis/growth/model_growth_censored.stan') # load growth model 

# STAN pars -------------- 
ncores <- 4 
niter <- 2000 
nchains <- 4 
nthin <- 4
# --------------------------
  
# pars used for all species ----------------
vr <- 'growth'
hold <- c(26:30)          ### Choose the hold out years 
formC <- as.formula(~-1)  ### Climate effects design matrix 
# ------------------------------------------

# Fit ARTR growth model -------------------------------------------- # 

dat_file <- 'data/temp_data/ARTR_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

dat$size <- scale( dat$logarea.t0 )
hist(dat$size) ## histogram of size to choose "small" size cutoff 

# pars ------------------- 
small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formE = as.formula(~ size)  ### For growth model, size dependent variance design matrix  
formX = as.formula(~ size + W + GroupP2 + Treatment2)  ### Fixed effects design matrix (include climate as "C")
# ---------------------- 

dat$small <- as.numeric(dat$size < small)
dat$Y    <- scale( dat$logarea.t1 )
dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 

dl <- process_data(dat = dat, 
                   formX = formX, 
                   formC = formC, 
                   formZ = formZ, 
                   formE = formE, 
                   center = T, 
                   vr = 'growth', 
                   hold )

# Assign a cut-off for data censoring of the Y-values 
freq <- data.frame( table(dl$Y) )
freq[ order(freq$Var1), ][1:10, ] 
hist(dl$Y)
left_cut <- -3 
dl <- left_censor(dl, U = left_cut)

# save the datalist 
saveRDS(dl, '~/Desktop/ARTR_growth_data.RDS')

gfit1 <- rstan::sampling(gmod, 
                         data = dl, 
                         chains = nchains, 
                         iter = niter, 
                         cores = ncores, 
                         thin = nthin )

saveRDS(gfit1, '~/Desktop/ARTR_growth_fit1.RDS')


## Fit HECO ------------------------------------------- # 

dat_file <- 'data/temp_data/HECO_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

dat$size <- scale( dat$logarea.t0 )
hist(dat$size) ## histogram of size to choose "small" size cutoff 

# pars ------------------- 
small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formE = as.formula(~ size)  ### For growth model, size dependent variance design matrix  
formX = as.formula(~ size + W + GroupP2 + Treatment2)  ### Fixed effects design matrix (include climate as "C")
# ---------------------- 

dat$small <- as.numeric(dat$size < small)
dat$Y    <- scale( dat$logarea.t1 )
dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 

dl <- process_data(dat = dat, 
                   formX = formX, 
                   formC = formC, 
                   formZ = formZ, 
                   formE = formE, 
                   center = T, 
                   vr = 'growth', 
                   hold )

# Assign a cut-off for data censoring of the Y-values 
freq <- data.frame( table(dl$Y) )
freq[ order(freq$Var1), ][1:15, ] 
hist(dl$Y)
left_cut <- -1.6 
dl <- left_censor(dl, U = left_cut)

# save the datalist 
saveRDS(dl, '~/Desktop/HECO_growth_data.RDS')

gfit1 <- rstan::sampling(gmod, 
                         data = dl, 
                         chains = nchains, 
                         iter = niter, 
                         cores = ncores, 
                         thin = nthin )

saveRDS(gfit1, '~/Desktop/HECO_growth_fit1.RDS')



## Fit POSE ------------------------------------------- # 

dat_file <- 'data/temp_data/POSE_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

dat$size <- scale( dat$logarea.t0 )
hist(dat$size) ## histogram of size to choose "small" size cutoff 

# pars ------------------- 
small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formE = as.formula(~ size)  ### For growth model, size dependent variance design matrix  
formX = as.formula(~ size + W + GroupP2 + Treatment2)  ### Fixed effects design matrix (include climate as "C")
# ---------------------- 

dat$small <- as.numeric(dat$size < small)
dat$Y    <- scale( dat$logarea.t1 )
dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 

dl <- process_data(dat = dat, 
                   formX = formX, 
                   formC = formC, 
                   formZ = formZ, 
                   formE = formE, 
                   center = T, 
                   vr = 'growth', 
                   hold )

# Assign a cut-off for data censoring of the Y-values 
freq <- data.frame( table(dl$Y) )
freq[ order(freq$Var1), ][1:15, ] 
hist(dl$Y)
left_cut <- -1.68 
dl <- left_censor(dl, U = left_cut)

# save the datalist 
saveRDS(dl, '~/Desktop/POSE_growth_data.RDS')

gfit1 <- rstan::sampling(gmod, 
                         data = dl, 
                         chains = nchains, 
                         iter = niter, 
                         cores = ncores, 
                         thin = nthin )

saveRDS(gfit1, '~/Desktop/POSE_growth_fit1.RDS')



## Fit PSSP ------------------------------------------- # 

dat_file <- 'data/temp_data/PSSP_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

dat$size <- scale( dat$logarea.t0 )
hist(dat$size) ## histogram of size to choose "small" size cutoff 

# pars ------------------- 
small <- -1               ### Designate a "small" size theshhold 
formZ = as.formula(~ size)  ### Year effects design matrix 
formE = as.formula(~ size)  ### For growth model, size dependent variance design matrix  
formX = as.formula(~ size + W + GroupP2 + Treatment2)  ### Fixed effects design matrix (include climate as "C")
# ---------------------- 

dat$small <- as.numeric(dat$size < small)
dat$Y    <- scale( dat$logarea.t1 )
dat$GroupP2 <- as.numeric( dat$Group == 'P2') # Paddock P2 is weird 

dl <- process_data(dat = dat, 
                   formX = formX, 
                   formC = formC, 
                   formZ = formZ, 
                   formE = formE, 
                   center = T, 
                   vr = 'growth', 
                   hold )

# Assign a cut-off for data censoring of the Y-values 
freq <- data.frame( table(dl$Y) )
freq[ order(freq$Var1), ][1:25, ] 
hist(dl$Y)
left_cut <- -1.95 
dl <- left_censor(dl, U = left_cut)

# save the datalist 
saveRDS(dl, '~/Desktop/PSSP_growth_data.RDS')

gfit1 <- rstan::sampling(gmod, 
                         data = dl, 
                         chains = nchains, 
                         iter = niter, 
                         cores = ncores, 
                         thin = nthin )

saveRDS(gfit1, '~/Desktop/PSSP_growth_fit1.RDS')

