##### Spring 2015 soil moisture
##### 

rm(list = ls())
library (ggplot2 )
library( reshape )

plots = read.csv('../precip_plots.csv')
plots = plots[ , c('plot_id', 'block', 'treatment') ]

readings = read.csv('2015_spring/2015-04-29-precip_experiment_soil_moisture.csv', skip = 2)
readings_long = melt(readings, id.var = 'plot')
moisture = merge( plots, readings_long, by.x = 'plot_id', by.y = 'plot' , all.y = TRUE)
levels(moisture$treatment ) <- c(levels(moisture$treatment), 'control')
moisture$treatment[ is.na(moisture$treatment)] <-  'control'

theme_set(theme_classic())
ggplot( moisture, aes( x = treatment, y = value )) + geom_boxplot() + 
  geom_point( position = position_jitter(width= 0.01)) + ylab('Percent Soil Moisture % VWC')

