rm(list = ls() )

library(rstan)
library(parallel)
library(dplyr)
library(tidyr)

# input ------------------------------------------------------------------------------------# 
lppd_table  <- read.csv('output/lppd_scores.csv')
lppd_table <- lppd_table %>% filter( model %in% c('treatment', 'climate', 'year'))

total_lppd <- 
  lppd_table %>% 
  group_by( species, vital_rate , model ) %>% 
  summarise( lppd = sum(lppd2))

overall_rank <- 
  total_lppd %>% 
  arrange( species, vital_rate, lppd )

rank_by_exp <- 
  lppd_table %>% 
  mutate( Experimental = ifelse (Treatment != 'Control', 'Experiment', 'Control')) %>% 
  group_by( species, vital_rate, Experimental, model ) %>% 
  summarise( lppd = sum(lppd2)) %>% 
  arrange( Experimental, species, vital_rate, lppd ) %>% 
  filter( lppd == min(lppd) )

rank_by_Treatment <- 
  lppd_table %>% 
  group_by( species, vital_rate, Treatment, model ) %>% 
  summarise( lppd = sum(lppd2)) %>% 
  arrange( vital_rate, species, Treatment, lppd ) %>%   
  filter( lppd == min(lppd) )

write.csv(overall_rank, 'output/model_scores.csv', row.names = FALSE)
#write.csv(best_models, 'output/lppd_selected_models.csv', row.names = FALSE)


