rm(list = ls())
library(rstan)

my_fit <- readRDS('output/stan_fits/ARTR_growth_treatment_fit.RDS')

treatment_effect_samples <- data.frame( rstan::extract(my_fit,'bt' ))

dat <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')


### Peter e-mail me these files as attachments and I will check that they match my output when 
### I run the model WITH the interaction (matching yours)

saveRDS( dat$ARTR, file = 'output/ARTR_growth_datalist.RDS')
saveRDS( treatment_effect_samples, file = 'output/ARTR_growth_treat_effects.RDS')


