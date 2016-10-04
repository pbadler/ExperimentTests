# test run on HPC  --------------------------------------------------------------------------------------------------- # 
rm(list = ls())

library(rstan)

setwd('/projects/A01633220/precip/')

datalist <- readRDS('data/temp_data/recruitment_data_lists_for_stan.RDS')
init_vals <- readRDS('data/temp_data/recruitment_init_vals.RDS')

test_spp <- 'ARTR'
test <- datalist[[test_spp]]
inits <- init_vals[[test_spp]]
model_no <- 5

#test$parents1 <- test$parents1 [ , grep( colnames (test$parents1), pattern = test_spp)  ] 
#test$parents2 <- test$parents2 [ , grep( colnames (test$parents2), pattern = test_spp)  ] 

#test$parents1_out <- test$parents1_out [ , grep( colnames (test$parents1_out), pattern = test_spp)  ] 
#test$parents2_out <- test$parents2_out [ , grep( colnames (test$parents2_out), pattern = test_spp)  ] 

inits <- rep(inits[model_no], 1)

test$tau_beta <- 1.5

test <- test[c('N', 'Nspp', 'spp',  'Y', 'gid', 'G', 'Yrs', 'yid', 'C', 'Covs', 'parents1', 'parents2', 'tau_beta', 
               'npreds', 'nyrs_out', 'yid_out', 'gid_out', 'y_holdout', 'Chold', 'parents1_out', 'parents2_out')]

init_fit <- readRDS('output/stan_fits/ARTR_recruitment_5_20_0_.RDS')

test_fit <- 
  stan(fit = init_fit, 
                 data = test, chains = 1, iter = 10, cores = 1)

test_fit <- stan(file = 'analysis/recruitment/model_recruitment_5.stan', data = test, chains = 1, iter = 10, cores = 1)


saveRDS(test_fit, 'output/stan_fits/predictions/test_recruitment_fit.RDS')
