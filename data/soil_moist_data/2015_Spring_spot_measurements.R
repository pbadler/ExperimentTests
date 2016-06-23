##### Spring 2015 soil moisture
##### 

rm(list = ls())
library(ggplot2 )
library(dplyr)
library(tidyr)

q_info <- read.csv('data/quad_info.csv')


calibration <- read.csv('data/soil_moist_data/spot_measurements/2015-05-07_soil_probe_calibration.csv')


p1 <- read.csv('data/soil_moist_data/spot_measurements/2012-06-06_spot_measurements.csv', skip = 3)

p1$date <- '2012-06-06'
p1$Plot <- gsub( p1$Plot, pattern = '-', replacement = '_')
p1$rep <- c(1:2)
p1 <- p1 %>% rename( plot = Plot )

p2 <- read.csv('data/soil_moist_data/spot_measurements/2015-05-07_spot_measurements.csv') 
p3 <- read.csv('data/soil_moist_data/spot_measurements/2016-05-10_spot_measurements.csv')
p4 <- read.csv('data/soil_moist_data/spot_measurements/2015-06-09_spot_measurements.csv')

df <- rbind( p2, p3, p4)

df <- df %>% gather( key = rep, PCT, E1:W3 )

df <- rbind( p1, df )

# run calibration ---------------------------------------------------------------------------------------------------- 

calibration %>% 
  mutate( wet_weight = (wet_weight- bag_weight) , dry_weight = (dry_weight - bag_weight), 
          water = wet_weight - dry_weight, soil_moisture_g = 100*water/dry_weight) 


readings$plot_id <- readings$plot

readings
merge( readings, plots, by= 'plot_id' ) 


readings_long = melt(readings, id.var = 'plot')
moisture = merge( plots, readings_long, by.x = 'plot_id', by.y = 'plot' , all.y = TRUE)
levels(moisture$treatment ) <- c(levels(moisture$treatment), 'control')
moisture$treatment[ is.na(moisture$treatment)] <-  'control'

theme_set(theme_classic())
ggplot( moisture, aes( x = treatment, y = value )) + geom_boxplot() + 
  geom_point( position = position_jitter(width= 0.01)) + ylab('Percent Soil Moisture % VWC')

