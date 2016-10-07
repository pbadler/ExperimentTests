rm(list = ls() )

library(rstan)
library(parallel)
library(dplyr)
library(tidyr)

# # function --------------------------------------------------------------------------------# 
# 
compute_lppd <- function( stan_fit ) {
  log_lik <- rstan::extract (stan_fit, "log_lik")$log_lik
  lppd <- log(colMeans(exp(log_lik)))
  sum(lppd)
}

# input ------------------------------------------------------------------------------------# 
waics <- read.csv('output/best_WAIC_scores.csv')
lppd  <- read.csv('output/lppd_scores.csv')

#real_files <- dir('output/stan_fits/predictions', pattern = '*[4]_predict*', full.names = TRUE)
#waics$pred_files <- file.path( 'output', 'stan_fits', 'predictions', paste( waics$species, waics$vital_rate, waics$model, waics$prior, waics$chains, 'predict.RDS', sep = '_'))
# waics$pred_files %in% real_files
# waics <- 
#   waics %>% 
#   filter( !vital_rate == 'recruitment' )
# lppd <- mclapply( waics$pred_files, FUN = function( x ) compute_lppd, mc.cores = 4)
# waics$prediction_lppd <- lppd

total_lppd <- 
  lppd %>% 
  filter( Period == 'Modern') %>% 
  group_by( species, vital_rate , model ) %>% 
  summarise( lppd_out = sum(lppd)) %>% 
  group_by( species, vital_rate) %>% 
  mutate( rank_lppd = lppd_out/sum(lppd_out))


model_scores <- merge( total_lppd, waics, by = c('species', 'vital_rate', 'model')  )

best_models <- 
  model_scores %>% 
  group_by(species, vital_rate) %>% 
  filter( lppd_out == max(lppd_out))

best_WAIC_models <- 
  waics %>% 
  group_by( species, vital_rate ) %>% 
  filter( waic == min(waic))

write.csv(best_WAIC_models, 'output/WAIC_selected_models.csv', row.names = FALSE)
write.csv(model_scores, 'output/model_scores.csv', row.names = FALSE)
write.csv(best_models, 'output/lppd_selected_models.csv', row.names = FALSE)
