rm(list = ls() )

library(ggmcmc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(rstan)
library(gridExtra)

my_colors <- c('#66c2a5','#fc8d62','#8da0cb')

# input ------------------------------------------------------------------------------------# 
model_table <- read.csv('output/WAIC_selected_models.csv')
lppd <- read.csv('output/lppd_scores.csv')

model_table <- model_table %>% filter( vital_rate != 'recruitment')

lppd_table <- merge(lppd, model_table, by = c('species', 'vital_rate' , 'model') )

avgs <- 
  lppd_table %>% 
  filter( Period == 'Modern') %>% 
  group_by(species, year, vital_rate, treatment) %>%
  summarise( obs_avg = mean(Y), pred_avg = mean(X50.))

ts_plot <- 
  ggplot(avgs, aes( x = year, y = obs_avg, color = treatment)) + 
  geom_point() + 
  geom_line() +
  scale_color_manual(values = my_colors) + 
  facet_wrap(~ species, nrow = 1)

ts_plot %+% filter( avgs , vital_rate == 'growth') 




