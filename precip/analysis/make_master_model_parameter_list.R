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
niter <- c(2000, 2000, 2000)

make_pars_string <- function( x ) { 
  paste0( 'c(', toString(sQuote( x )), ')')
}

survparms1 <- c('a','b1_mu','b1','w','b2', 'mu','muhat','log_lik','log_lik2', 'bg','sig_a','sig_b1')
growparms1 <- c('a','b1_mu','b1','w','b2', 'mu','muhat','log_lik','log_lik2','bg','sig_a','sig_b1','sigma', 'sigmahat', 'tau', 'tauSize')
recparms1 <- c('a','theta','u','w', 'b2','log_lik','log_lik2','lambda', 'lambda_pred', 'bg', 'sig_a')

smodels <- data.frame(model = 1, 
                      vital_rate = 'survival',
                      pars=c( make_pars_string(survparms1))) 

gmodels <-data.frame( model = 1, 
                      vital_rate = 'growth',  
                      pars = c(make_pars_string(growparms1)))

rmodels <- data.frame(model = 1, 
                      vital_rate = 'recruitment', 
                      pars = c( make_pars_string(recparms1)))


# regularization based on Gerber et al. 2015 ---------------------------------------------------------------------# 
nlambda <- 40
print( paste ( 'nlambda =' , nlambda ) ) 
lambda.set <- exp(seq(-4.5, 10, length.out = nlambda))
sd_vec <- sqrt(1/lambda.set) # use sd for stan normal distribution 

sd_vec <- cbind( lambda = 1:nlambda, sd = sd_vec )

# ----------------------------------------------------------------------------------------------------------------# 

ff <- dir('data/temp_data', '*_*_cleaned_dataframe.RDS', full.names = T)

yrs <- lapply(ff,  function(x) unique( readRDS(x)$year )) 
sps <- str_extract(basename(ff), '^[A-Z]+')
vrs <- str_extract(basename(ff), '(growth)|(recruitment)|(survival)')

dfyrs <- list(NA)
for( i in 1:length(sps)) { 
   dfyrs[[i]] <- data.frame( species = sps[i], vital_rate = vrs[i], year_oos = yrs[[i]] )
}
dfyrs <- do.call(rbind, dfyrs)

dfyrs$year_oos[ dfyrs$year_oos > 2006 ] <- 'c(2007:2015)'

dfyrs <- unique(dfyrs)

vital_rates <- data.frame( vital_rate = vital_rates, niter = niter)

master_list <- list(species = species, vital_rates = vital_rates, smodels = smodels, gmodels = gmodels, rmodels = rmodels, nlambda = nlambda, sd_vec = sd_vec, dfyrs = dfyrs)

save(master_list, file = 'data/temp_data/master_list.Rdata')

