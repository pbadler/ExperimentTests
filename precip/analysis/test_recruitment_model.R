library(rstan)
models <- read.csv('data/temp_data/model_table.csv')

datalist <- readRDS('data/temp_data/recruitment_data_lists_for_stan.RDS')[['ARTR']]

init_vals <- readRDS( 'data/temp_data/recruitment_init_vals.RDS')[['ARTR']]

init_vals[3]

sum( datalist$parents1 < 0 ) 
sum( datalist$parents2 < 0 )

fit <- stan(file = 'analysis/recruitment/model_recruitment_3_test.stan', data= datalist, chains = 1, iter = 1)

fit <- stan(file = 'analysis/recruitment/test_stan_sqrt.stan', data = datalist, chains = 1, iter = 1)
