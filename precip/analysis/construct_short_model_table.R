# 
#
# construct initialize model table 
#
#
library(dplyr)
library(tidyr)

rm(list = ls())

do_spp <- c('ARTR','HECO', 'POSE', 'PSSP')
do_vr_number <- c(1, 2, 3)
do_vr <-  c('growth', 'survival', 'recruitment')
do_model <- c(1:5)
prior <- 15

models <- expand.grid(species = do_spp, vital_rate_number = do_vr_number, model = do_model, prior  = prior ) 

vital_rates <- data.frame(vital_rate = do_vr, vital_rate_number = do_vr_number)

model_type <- data.frame( model = do_model, pars = c( "c('a_mu', 'b1_mu', 'mu', 'muhat',  'log_lik', 'y_hat')", 
                                                 "c('a_mu', 'b1_mu', 'b2', 'mu', 'muhat', 'log_lik', 'y_hat')", 
                                                 "c('a_mu', 'b1_mu', 'w', 'mu', 'muhat',  'log_lik', 'y_hat')", 
                                                 "c('a_mu', 'b1_mu', 'b2', 'w', 'mu', 'muhat',  'log_lik', 'y_hat')", 
                                                 "c('a_mu', 'b1_mu', 'b2', 'w', 'mu', 'muhat',  'log_lik', 'y_hat')"))

models <- merge( models, model_type)

models <- merge(models, vital_rates)

models <- 
  models %>% 
  arrange( vital_rate_number, species, model) %>% 
  dplyr::select( vital_rate_number, vital_rate, species, model, prior , pars )

# --------------------------------------------------------------------------------------- 

write.csv( models, 'data/temp_data/short_model_table.csv', row.names = FALSE)

# ---------------------------------------------------------------------------------------

