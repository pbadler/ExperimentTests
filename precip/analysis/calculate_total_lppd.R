rm(list = ls() )

library(rstan)
library(parallel)
library(dplyr)
library(tidyr)

# input ------------------------------------------------------------------------------------# 
waics <- read.csv('output/best_WAIC_scores.csv')
lppd_table  <- read.csv('output/lppd_scores.csv')

total_lppd <- 
  lppd_table %>% 
  group_by( species, vital_rate , model ) %>% 
  summarise( lppd_out = sum(lppd1), lppd_out2 = sum(lppd2)) %>% 
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
