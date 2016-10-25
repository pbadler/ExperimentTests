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
niter <- c(4000, 10000, 4000)

smodels <- data.frame(model = 1:3, 
                      vital_rate = 'survival',
                      pars=c(  "c('a',  'b1_mu', 'b1',  'w',  'mu',  'muhat',  'log_lik', 'y_hat', 
                                  'a2', 'b1_mu2', 'b12', 'w2', 'mu2', 'muhat2', 'log_lik2','y_hat2')", 
                               "c('a',  'b1_mu', 'b12',  'b2', 'w',   'mu',     'muhat',   'log_lik', 'y_hat', 
                                  'a2', 'b1_mu2', 'b12', 'w2', 'mu2', 'muhat2', 'log_lik2','y_hat2')", 
                               "c('a',  'b1_mu', 'b12', 'b2', 'w',   'mu',     'muhat',   'log_lik', 'y_hat',
                                  'a2', 'b1_mu2', 'b12', 'b2', 'w2', 'mu2', 'muhat2', 'log_lik2','y_hat2')"))

gmodels <-data.frame( model = 1:3, 
                      vital_rate = 'growth',  
                      pars = c("c('a_mu',  'a',  'b1_mu', 'b1', 'w', 'mu', 'muhat',  'log_lik', 'y_hat', 
                                  'a_mu2', 'a2', 'b1_mu2', 'b12', 'w2', 'mu2', 'muhat2',  'log_lik2', 'y_hat2', 
                                  'muhat3','y_hat3', 
                                  'muhat4','y_hat4')", 
                               "c('a_mu',  'a',  'b1_mu', 'b1', 'b2', 'w', 'mu', 'muhat',  'log_lik', 'y_hat', 
                                  'a_mu2', 'a2', 'b1_mu2', 'b2', 'w2', 'mu2', 'muhat2',  'log_lik2', 'y_hat2', 
                                  'muhat3','y_hat3', 
                                  'muhat4','y_hat4')", 
                               "c('a_mu',  'a',  'b1_mu', 'b1', 'b2', 'w', 'mu', 'muhat',  'log_lik', 'y_hat', 
                                  'a_mu2', 'a2', 'b1_mu2', 'b12', 'w2', 'mu2', 'muhat2',  'log_lik2', 'y_hat2', 
                                  'muhat3','y_hat3', 
                                  'muhat4','y_hat4')"))

rmodels <- data.frame(model = 1:3, 
                      vital_rate = 'recruitment', 
                      pars = c( "c('a_mu', 'a',    'theta',  'u',  'q',  'w',  'qpred',  'log_lik',  'y_hat', 
                                  'a_mu2', 'a2',   'theta2', 'u2', 'q2', 'w2', 'qpred2', 'log_lik2', 'y_hat2')", 
                                "c('a_mu', 'a',    'theta',  'u',  'q',  'b2', 'w',      'qpred',    'log_lik', 'y_hat', 
                                  'a_mu2', 'a2',   'theta2', 'u2', 'q2', 'w2', 'qpred2', 'log_lik2', 'y_hat2')", 
                                "c('a_mu', 'a',    'theta',  'u',  'q',  'b2', 'w',      'qpred',    'log_lik', 'y_hat', 
                                  'a_mu2', 'a2',   'theta2', 'u2', 'q2', 'w2', 'qpred2', 'log_lik2', 'y_hat2')"))


# regularization based on Gerber et al. 2015 ---------------------------------------------------------------------# 
nlambda <- 40
print( paste ( 'nlambda =' , nlambda ) ) 
lambda.set <- exp(seq(-4.5, 10, length.out = nlambda))
sd_vec <- sqrt(1/lambda.set) # use sd for stan normal distribution 

sd_vec <- cbind( lambda = 1:nlambda, sd = sd_vec )

# ----------------------------------------------------------------------------------------------------------------# 

vital_rates <- data.frame( vital_rate = vital_rates, niter = niter)

master_list <- list(species = species, vital_rates = vital_rates, smodels = smodels, gmodels = gmodels, rmodels = rmodels, nlambda = nlambda, sd_vec = sd_vec)

save(master_list, file = 'data/temp_data/master_list.Rdata')

