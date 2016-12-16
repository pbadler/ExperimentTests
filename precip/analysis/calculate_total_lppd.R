rm(list = ls() )

library(rstan)
library(parallel)
library(dplyr)
library(tidyr)

# input ------------------------------------------------------------------------------------# 
lppd_table  <- read.csv('output/lppd_scores.csv')
lppd_table <- lppd_table %>% filter( model %in% c('treatment', 'climate', 'year'))

overall_rank <- 
  lppd_table %>% 
  group_by( species, vital_rate , model ) %>% 
  summarise( lppd = sum(lppd2)) %>% 
  arrange( species, vital_rate, lppd )

rank_by_exp <- 
  lppd_table %>% 
  mutate( Experimental = ifelse (Treatment != 'Control', 'Experiment', 'Control')) %>% 
  group_by( species, vital_rate, Experimental, model ) %>% 
  summarise( lppd = sum(lppd2)) %>% 
  arrange( Experimental, species, vital_rate, lppd ) %>% 
  group_by( Experimental, species, vital_rate ) %>% 
  spread( model, lppd ) %>% 
  mutate( climate_model_performance = climate - year )

rank_by_Treatment <- 
  lppd_table %>% 
  group_by( species, vital_rate, Treatment, model ) %>% 
  summarise( lppd = sum(lppd2)) %>% 
  arrange( vital_rate, species, Treatment, lppd ) %>%  
  group_by( Treatment, species, vital_rate ) %>% 
  spread( model, lppd ) %>% 
  mutate( climate_model_performance = climate - year )

rank_by_size_class <- 
  lppd_table %>%
  filter( vital_rate != 'recruitment' ) %>% 
  group_by( species, vital_rate, model ) %>% 
  mutate( size_class  = cut( size , 3 , labels = c('s', 'm' ,'l')) ) %>% 
  group_by( species, vital_rate, model, size_class ) %>% 
  summarise( lppd  = sum(lppd2 ) ) %>% 
  arrange( size_class, species, vital_rate, lppd ) %>% 
  group_by( species, vital_rate, size_class) %>% 
  spread( model, lppd ) %>% 
  mutate( climate_model_performance = climate - year )

rank_by_Treatment

overall_rank


write.csv(overall_rank, 'output/model_scores.csv', row.names = FALSE)
#write.csv(best_models, 'output/lppd_selected_models.csv', row.names = FALSE)


