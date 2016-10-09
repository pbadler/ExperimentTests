rm(list = ls())
library(rstan)

models <- read.csv(file = 'data/temp_data/model_table.csv')

datalist <- readRDS('data/temp_data/recruitment_data_lists_for_stan.RDS')[['PSSP']]

datalist_old <- readRDS('data/temp_data/recruitment_data_that_works.RDS')[['PSSP']]

datalist$Y
datalist_old$Y

head( datalist$parents1 ) 
head( datalist_old$parents1)

identical(datalist_old$parents1, datalist$parents1)
identical(datalist_old$parents2, datalist$parents2)

head( datalist$parents2 ) 
head( datalist_old$parents2)

datalist$tau_beta <- 1000
sqrt( 1/(exp(10)))
sqrt( 1/(exp(-5)) ) 

init_vals <- readRDS( 'data/temp_data/recruitment_init_vals.RDS')[['PSSP']]
init_vals <- rep(init_vals[3], 4)
fit <- stan(file = 'analysis/recruitment/model_recruitment_3.stan', init = init_vals, data= datalist, chains = 4, iter = 2000, cores = 4)

summary(fit, c('a_mu','sig_a', 'sig_G', 'a' , 'gint', 'theta', 'u', 'w', 'b2'))$summary

traceplot(fit, 'b2')

datalist$Y[242]
datalist$parents1[242, ]
datalist$parents2[242, ]
datalist$C[273, ]
matplot(datalist$C[200:300, ])
