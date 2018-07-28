rm(list = ls())

library(tidyverse)

grid <- expand.grid(species = c('ARTR', 'HECO', 'POSE', 'PSSP'), vr = c( 'survival', 'growth'))


for( i in 1:nrow(grid)){ 
  spp <- grid$sp[i]
  vr  <- grid$vr[i]
  scores <- readRDS(paste0( 'output/', spp, '_', vr, '_model_scores.RDS'))
  
  scores %>% 
    ggplot(aes(x = climate_effects, y = out_of_sample_lppd)) + 
    geom_point() 
  
  scores %>% 
    ggplot(aes(x = climate_effects, y = out_of_sample_mse)) + 
    geom_point() 
  

  write_csv( scores %>% arrange(desc(out_of_sample_lppd)), paste0('output/', spp, '_', vr, '_', 'model_ranks.csv'))

}

