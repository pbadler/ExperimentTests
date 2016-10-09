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
niter <- c(2000, 5000, 2000)

models <- data.frame( model = 1:3, pars = c("c('a_mu', 'b1_mu', 'w', 'mu', 'muhat',  'log_lik', 'y_hat')", 
                                              "c('a_mu', 'b1_mu', 'b2', 'w', 'mu', 'muhat',  'log_lik', 'y_hat')", 
                                              "c('a_mu', 'b1_mu', 'b2', 'w', 'mu', 'muhat',  'log_lik', 'y_hat')"))

rmodels <- data.frame( model = 1:3, pars = c("c('a_mu', 'theta', 'u', 'q', 'w', 'qpred', 'log_lik', 'y_hat')", 
                                             "c('a_mu', 'theta', 'u', 'q', 'b2', 'w', 'qpred',  'log_lik', 'y_hat')", 
                                             "c('a_mu',  'theta', 'u', 'q', 'b2', 'w', 'qpred',  'log_lik', 'y_hat')"))

# regularization based on Gerber et al. 2015 ---------------------------------------------------------------------# 
nlambda <- 40
print( paste ( 'nlambda =' , nlambda ) ) 
lambda.set <- exp(seq(-5, 10, length=nlambda))
sd_vec <- sqrt(1/lambda.set) # use sd for stan normal distribution 

sd_vec <- cbind( lambda = 1:nlambda, sd = sd_vec )

# ----------------------------------------------------------------------------------------------------------------# 

vital_rates <- data.frame( vital_rate = vital_rates, niter = niter)

master_list <- list(species = species, vital_rates = vital_rates, models = models, rmodels = rmodels, nlambda = nlambda, sd_vec = sd_vec)

save(master_list, file = 'data/temp_data/master_list.Rdata')

