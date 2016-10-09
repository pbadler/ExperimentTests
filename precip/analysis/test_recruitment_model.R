library(rstan)
models <- read.csv('data/temp_data/model_table.csv')

datalist <- readRDS('data/temp_data/recruitment_data_lists_for_stan.RDS')[['POSE']]

init_vals <- readRDS( 'data/temp_data/recruitment_init_vals.RDS')[['POSE']]

fit <- stan(file = 'analysis/recruitment/model_recruitment_3.stan', data= datalist, init = rep(init_vals[3], 4), chains = 4, iter = 10000)


summary(fit, c('a_mu','sig_a', 'sig_G', 'a' , 'gint', 'theta', 'u', 'w', 'b2'))
