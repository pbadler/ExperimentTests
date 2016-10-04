# #######################################################################################
#
# Construct master model parameter list 
# list includes 
# all species
# all vital rates 
# all models 
# number of priors 
#
# #######################################################################################

rm(list = ls())

library(tidyr)
library(dplyr)

species <- c('ARTR', 'HECO', 'POSE', 'PSSP')
vital_rates <- c('growth', 'recruitment', 'survival')

models <- data.frame( model = 1:5, pars = c( "c('a_mu', 'b1_mu', 'mu', 'muhat',  'log_lik', 'y_hat')", 
                                              "c('a_mu', 'b1_mu', 'b2', 'mu', 'muhat', 'log_lik', 'y_hat')", 
                                              "c('a_mu', 'b1_mu', 'w', 'mu', 'muhat',  'log_lik', 'y_hat')", 
                                              "c('a_mu', 'b1_mu', 'b2', 'w', 'mu', 'muhat',  'log_lik', 'y_hat')", 
                                              "c('a_mu', 'b1_mu', 'b2', 'w', 'mu', 'muhat',  'log_lik', 'y_hat')"))

rmodels <- data.frame( model = 1:5, pars = c("c('a_mu', 'theta', 'u', 'q', 'log_lik', 'qpred', 'y_hat')", 
                                             "c('a_mu', 'theta', 'u', 'q', 'b2', 'qpred', 'log_lik', 'y_hat')", 
                                             "c('a_mu', 'theta', 'u', 'q', 'w', 'qpred', 'log_lik', 'y_hat')", 
                                             "c('a_mu', 'theta', 'u', 'q', 'b2', 'w', 'qpred',  'log_lik', 'y_hat')", 
                                             "c('a_mu', 'b1_mu', 'u', 'q', 'b2', 'w', 'qpred',  'log_lik', 'y_hat')"))

nlambda <- 40

master_list <- list(species = species, vital_rates = vital_rates, models = models, rmodels = rmodels, nlambda = nlambda)

save(master_list, file = 'data/temp_data/master_list.Rdata')