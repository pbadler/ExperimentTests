rm(list = ls() )

library(rstan)
library(parallel)
library(dplyr)

# function --------------------------------------------------------------------------------# 

compute_lppd <- function( stan_fit ) { 
  log_lik <- extract (stan_fit, "log_lik")$log_lik
  lppd <- log(colMeans(exp(log_lik)))
  sum(lppd)
} 

# input ------------------------------------------------------------------------------------# 
setwd('~/Documents/ExperimentTests/precip/')
print(paste('Working directory: ' , getwd()))

waics <- readRDS('output/best_WAIC_scores.RDS')

waics <- waics %>% 
  mutate( prediction_file = paste( species, vital_rate, model, 'predictions.RDS', sep = '_' )) 

pred_files <- dir('output/stan_fits/predictions/', 'predictions.RDS', full.names = TRUE)

lppd <- mclapply( pred_files , function(x ) compute_lppd( readRDS(x) ), mc.cores = 4)

prediction_results <- data.frame(prediction_file = basename(pred_files), prediction_lppd = unlist ( lppd)  )

prediction_results <- merge( waics, prediction_results, by = 'prediction_file', all.x = TRUE)

saveRDS(prediction_results, 'output/prediction_lppds.RDS')
