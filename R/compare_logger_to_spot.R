rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(lme4)
library(zoo)

df <- readRDS('data/temp_data/decagon_data_with_rainfall_data.RDS')
df_spot <- readRDS('data/temp_data/spring_spot_measurements.RDS')

# ----------------------------------------------------------------------------- 

head(df)

head( df_spot ) 

class(df$date)
class(df_spot$date)


daily_VWC <- df %>% 
  filter ( measure == 'VWC', depth == '5 cm deep') %>% 
  group_by ( plot, date, depth  ) %>% 
  summarise( logger_avg = mean( v ))

spot_avg <- 
  df_spot %>% 
  group_by( plot, date, Treatment, PrecipGroup) %>% 
  summarise( spot = mean(VWC))

daily_VWC$date <- as.Date(daily_VWC$date, tz = 'MST')
spot_avg$date <- as.Date( spot_avg$date, tz = 'MST')

spot_avg <- left_join( spot_avg, daily_VWC, by = c('plot', 'date'))

m1 <- lm(data = spot_avg, logger_avg ~ spot )
summary(m1)

#spot_avg <- subset( spot_avg, date > as.Date('2013-01-01') ) 

ggplot ( data = spot_avg, aes ( x = spot, y = logger_avg, color = factor(date) ) ) + geom_point()
         