rm(list = ls() )

library(ggmcmc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(rstan)
library(gridExtra)

# input ------------------------------------------------------------------------------------# 
model_table <- read.csv('output/WAIC_selected_models.csv')
lppd <- read.csv('output/lppd_scores.csv')

lppd_table <- merge(lppd, model_table, by = c('species', 'vital_rate' , 'model') )

lppd_table <- 
  lppd_table %>% 
  group_by( vital_rate, species, model ) %>% 
  mutate( diff = X50. - Y)

lppd_table %>% 
  group_by( vital_rate, species, model ) %>% 
  summarise( n())

g1 <-
  ggplot(lppd_table, aes(x = diff, fill = Period)) +
  geom_density(alpha = 0.6, color = NA) +
  geom_vline( aes(xintercept = 0), linetype = 2, alpha = 0.8) +
  xlab( 'Predicted - observed') + 
  facet_wrap( ~ Period, ncol = 1) 

plot1 <- lppd_table %>% do(gg = g1 %+% . + ggtitle( paste( .$species, .$vital_rate, 'model', .$model)))

my_colors <- c('#66c2a5','#fc8d62','#8da0cb')

gg_treatments <-
  ggplot( lppd_table, aes( x = diff, fill = treatment)) +
  geom_density(alpha = 0.6, color = NA) +
  facet_wrap(~ treatment, ncol = 1 ) +
  geom_vline( aes(xintercept = 0), linetype = 2, alpha = 0.8) +
  scale_fill_manual(values = my_colors) +
  xlab( 'Predicted - observed') 

plot2 <- 
  lppd_table %>% 
  filter( Period == 'Modern') %>% 
  do(gg = gg_treatments %+% . + ggtitle( paste( .$species, .$vital_rate, 'model', .$model)))

gg_years <-
  ggplot( lppd_table, aes( x = diff)) +
  geom_density(fill = 'gray', alpha = 0.6) +
  facet_wrap( ~ year, ncol = 1 ) +
  geom_vline( aes( xintercept = 0 ) , linetype = 2, alpha = 0.5 )

plot3 <- 
  lppd_table %>% 
  filter( Period == 'Modern') %>% 
  do(gg = gg_years %+% . + ggtitle( paste( .$species, .$vital_rate, 'model', .$model)))


png(file.path( 'figures/predictions/'))

for( i in 1:length(plot1$gg) ) { 
  
  vr <- plot1$vital_rate[i]
  species <- plot1$species[i]
  m <- plot1$model[i]
  
  png( file.path( 'figures/predictions', paste(species, vr, m, 'predicted_minus_obs.png', sep = '_')), width = 6, height = 6 , units = 'in', res = 300)
    print( plot1$gg[[i]] )
  dev.off()
  
  png( file.path( 'figures/predictions', paste(species, vr, m, 'predicted_minus_obs_by_treatment.png', sep = '_')), width = 6, height = 6 , units = 'in', res = 300)
    print( plot2$gg[[i]] )
  dev.off()
  
  png( file.path( 'figures/predictions', paste(species, vr, m, 'predicted_minus_obs_by_year.png', sep = '_')), width = 6, height = 6 , units = 'in', res = 300)
    print( plot3$gg[[i]] )
  dev.off()
  
  
}

dev.off()
