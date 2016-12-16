rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

files <- dir( 'output', 'treatment_model_parameters.*.csv', full.names = T)
load('analysis/figure_scripts/my_plotting_theme.Rdata')

dat <- lapply( files, read.csv)

dat <- do.call(rbind, dat )

treat <- dat[ str_detect(dat$X, 'bt'),  ] 

treat$significant <- ifelse(  ( treat$X2.5. > 0 & treat$X97.5. > 0 ) | (treat$X2.5. < 0 & treat$X97.5. < 0 ), '*', ' ' )

par_names <- data.frame( par = c('bt[1]', 'bt[2]', 'bt[3]', 'bt[4]'), par_name = c('Drought', 'Irrigation', 'Droughtxlogarea.t0', 'Irrigationxlogarea.t0'))
treat$par <- treat$X
treat <- merge( treat, par_names)
treat <- treat[order(treat$species, treat$vital_rate, treat$par_name), ]

treat$lbcl95 <- treat$X2.5. 
treat$ubcl95 <- treat$X97.5.
treat$Treatment <- NA
treat$Treatment[grep( 'D', treat$par_name )] <- 'Drought'
treat$Treatment[grep( 'I', treat$par_name )] <- 'Irrigation'
treat$type <- NA
treat$type[ grep( 'x', treat$par_name ) ] <- 'Treatment x size'
treat$type[ -grep( 'x', treat$par_name ) ] <- 'Treatment intercept'

treat <- treat[, c('species', 'vital_rate', 'Treatment', 'type', 'par_name', 'mean', 'se_mean', 'lbcl95', 'ubcl95', 'significant')]

effect_plot <- 
  ggplot( subset( treat, vital_rate == 'growth'), aes( x = Treatment, y = mean , ymin = lbcl95, ymax = ubcl95, color = Treatment )) + 
  geom_point() + 
  geom_errorbar() + 
  geom_hline( aes( yintercept = 0 ), linetype = 2 , alpha = 0.5) + 
  facet_grid( type  ~  species ) + 
  ylab ( 'mean effect (+/- 95% Bayesian Credible Interval)') + 
  scale_color_manual(values = my_colors[3:4]) + 
  my_theme

gg <- treat %>% group_by(vital_rate) %>% do(gg =  effect_plot %+% .  + ggtitle( paste0('Effects on ', .$vital_rate)))

pdf( 'figures/treatment_effect_estimates.pdf', height = 5, width = 8 ) 
print( gg$gg ) 
dev.off() 

