rm(list = ls())
library(rstan)

source('analysis/stan_data_functions.R')

mod <- rstan::stan_model('analysis/recruitment/recruitment.stan') # load survival model 

# STAN pars -------------- 
ncores <- 4 
niter <- 2000 
nchains <- 4 
nthin <- 4
# --------------------------

# pars used for all species ----------------
vr <- 'recruitment'
hold <- c(26:30)          ### Choose the hold out years 
formC <- as.formula(~-1)  ### Climate effects design matrix 
# ------------------------------------------

# Fit ARTR recruitment model -------------------------------------------- # 

dat_file <- 'data/temp_data/ARTR_recruitment_dataframe.RDS'
dat <- readRDS(dat_file)

# pars ------------------- 
formZ = as.formula(~ 1)  ### Year effects design matrix 
formX = as.formula(~ Treatment2)  ### Fixed effects design matrix (include climate as "C")
# ---------------------- 

dl <- process_recruitment_data(dat = dat, 
                   formX = formX, 
                   formC = formC, 
                   formZ = formZ, 
                   center = T, 
                   hold = hold )

# save the datalist 
saveRDS(dl, '~/Desktop/ARTR_recruitment_data.RDS')

rfit1 <- rstan::sampling(mod, 
                         data = dl, 
                         chains = nchains, 
                         iter = niter, 
                         cores = ncores, 
                         thin = nthin )


saveRDS(rfit1, '~/Desktop/ARTR_recruitment_fit1.RDS')


# Fit HECO recruitment model -------------------------------------------- # 

dat_file <- 'data/temp_data/HECO_recruitment_dataframe.RDS'
dat <- readRDS(dat_file)

# pars ------------------- 
formZ = as.formula(~ 1)  ### Year effects design matrix 
formX = as.formula(~ Treatment2)  ### Fixed effects design matrix (include climate as "C")
# ---------------------- 

dl <- process_recruitment_data(dat = dat, 
                               formX = formX, 
                               formC = formC, 
                               formZ = formZ, 
                               center = T, 
                               hold = hold )

# save the datalist 
saveRDS(dl, '~/Desktop/HECO_recruitment_data.RDS')

rfit1 <- rstan::sampling(mod, 
                         data = dl, 
                         chains = nchains, 
                         iter = niter, 
                         cores = ncores, 
                         thin = nthin )


saveRDS(rfit1, '~/Desktop/HECO_recruitment_fit1.RDS')



# Fit POSE recruitment model -------------------------------------------- # 

dat_file <- 'data/temp_data/POSE_recruitment_dataframe.RDS'
dat <- readRDS(dat_file)

# pars ------------------- 
formZ = as.formula(~ 1)  ### Year effects design matrix 
formX = as.formula(~ Treatment2)  ### Fixed effects design matrix (include climate as "C")
# ---------------------- 

dl <- process_recruitment_data(dat = dat, 
                               formX = formX, 
                               formC = formC, 
                               formZ = formZ, 
                               center = T, 
                               hold = hold )

# save the datalist 
saveRDS(dl, '~/Desktop/POSE_recruitment_data.RDS')

rfit1 <- rstan::sampling(mod, 
                         data = dl, 
                         chains = nchains, 
                         iter = niter, 
                         cores = ncores, 
                         thin = nthin )


saveRDS(rfit1, '~/Desktop/POSE_recruitment_fit1.RDS')



# Fit PSSP recruitment model -------------------------------------------- # 

dat_file <- 'data/temp_data/PSSP_recruitment_dataframe.RDS'
dat <- readRDS(dat_file)

# pars ------------------- 
formZ = as.formula(~ 1)  ### Year effects design matrix 
formX = as.formula(~ Treatment2)  ### Fixed effects design matrix (include climate as "C")
# ---------------------- 

dl <- process_recruitment_data(dat = dat, 
                               formX = formX, 
                               formC = formC, 
                               formZ = formZ, 
                               center = T, 
                               hold = hold )

# save the datalist 
saveRDS(dl, '~/Desktop/PSSP_recruitment_data.RDS')

rfit1 <- rstan::sampling(mod, 
                         data = dl, 
                         chains = nchains, 
                         iter = niter, 
                         cores = ncores, 
                         thin = nthin )


saveRDS(rfit1, '~/Desktop/PSSP_recruitment_fit1.RDS')


