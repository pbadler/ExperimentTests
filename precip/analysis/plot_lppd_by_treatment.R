rm(list = ls())

library(ggplot2 )
library(tidyr)
library(dplyr)

best_models <- read.csv('output/best_predicting_models.csv')

lppd <- read.csv('output/lppd_scores.csv')

best_lppd <- merge( lppd, best_models[, c('species', 'vital_rate', 'model')], by = c('species', 'vital_rate', 'model'))

lppd_by_treatment <- 
  best_lppd %>% 
  group_by( species, vital_rate, model, treatment) %>% 
  summarise( lppd =  mean(lppd)) %>%
  group_by( treatment) %>% 
  mutate( mlppd = mean(lppd))

png( 'figures/plot_lppd_by_treatment.png', 6, 6, units = 'in', res = 300 )
  print( ggplot( lppd_by_treatment, aes( x = treatment, y = lppd , fill = treatment) ) + 
    geom_boxplot() )
dev.off()

