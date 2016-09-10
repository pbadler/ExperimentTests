##########################################################################
#
# import datalist and fit STAN model 
# output stan model  
# 
##########################################################################

rm(list = ls())

library(rstan)

data_lists <- readRDS('data/temp_data/data_lists_for_stan_models.RDS')

pars <- c("log_lik")

stan_fit <- stan('analysis/growth/growth_oos_cv.stan', 'growth', data = data_lists[[1]], chains = 0, pars = pars )

