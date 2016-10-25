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
options("useFancyQuotes" = FALSE)
library(tidyr)
library(dplyr)

species <- c('ARTR', 'HECO', 'POSE', 'PSSP')
vital_rates <- c('growth', 'recruitment', 'survival')
niter <- c(4000, 10000, 4000)

make_pars_string <- function( x ) { 
  paste0( 'c(', toString(sQuote( x )), ')')
}

survparms1 <- c('a','b1_mu','b1','w','mu','muhat','log_lik','y_hat','bg','sig_a','sig_b1')
survparms1 <- c(survparms1, paste0(survparms1, '2'))
survparms2 <- c(survparms1, 'b2')

growparms1 <- c('a','b1_mu','b1','w','mu','muhat','log_lik','y_hat','bg','sig_a','sig_b1','sigma')
growparms1 <- c(growparms1, paste0(growparms1, '2'))
growparms1 <- c(growparms1, 'muhat3', 'y_hat3', 'muhat4', 'y_hat4')
growparms2 <- c(growparms1, 'b2')
growparms2 <- c(growparms2, 'muhat3', 'y_hat3', 'muhat4', 'y_hat4')

recparms1 <- c('a','theta','u','w','log_lik','q_pred', 'bg', 'sig_a')
recparms1 <- c(recparms1, paste0(recparms1, '2'))
recparms2 <- c(recparms1, 'b2')

smodels <- data.frame(model = 1:3, 
                      vital_rate = 'survival',
                      pars=c( make_pars_string(survparms1), make_pars_string(survparms2), make_pars_string(survparms2) )) 

gmodels <-data.frame( model = 1:3, 
                      vital_rate = 'growth',  
                      pars = c(make_pars_string(growparms1), make_pars_string(growparms2), make_pars_string(growparms2)))


rmodels <- data.frame(model = 1:3, 
                      vital_rate = 'recruitment', 
                      pars = c( make_pars_string(recparms1), make_pars_string(recparms2), make_pars_string(recparms2)))


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

