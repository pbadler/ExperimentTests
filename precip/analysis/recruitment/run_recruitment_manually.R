rm(list = ls())

load('data/temp_data/master_list.Rdata')
scores <- read.csv('output/best_WAIC_scores.csv')

spp <- 'HECO'
m <- 2
df <- readRDS('data/temp_data/recruitment_data_lists_for_stan.RDS')[[spp]]

library(rstan)

lambda <- subset(scores, species == spp & vital_rate == 'recruitment' & model == m)$lambda
tau_beta <- subset(scores, species == spp & vital_rate == 'recruitment' & model == m)$sd
df$tau_beta <- tau_beta

inits <- rep(list(list(w = rep(0, 1), b2 = rep(0, 7))), 4)

testfit1 <- stan(paste0( 'analysis/recruitment/model_recruitment_', m , '_predict.stan'), chains = 4, init = inits, data = df, cores = 4, iter = 2000)

# ----------------------------------------------------------------------------------------- # 

saveRDS(testfit1, paste0( 'output/stan_fits/predictions/', spp, '_recruitment_', m, '_', lambda, '_', '4_predict.RDS'))

# ----------------------------------------------------------------------------------------- # 