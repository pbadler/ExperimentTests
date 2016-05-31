rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)

df <- readRDS(file = 'data/temp_data/decagon_data.RDS')

# Correct bad dates ------------------------------------------------------------------------------ 

df %>% 
  filter( date < strptime('2012-01-01', '%Y-%m-%d') | date > strptime('2016-06-01', '%Y-%m-%d')) %>% 
  distinct(id, plot , date_started)

length( df$date[df$plot == 11 & df$date_started == '2014-04-13']  ) 

length( df$date[df$plot == '11_12_C' & df$date_started == '2014-04-13']  ) 

times <- df %>% select( PrecipGroup, period, plot, date ) %>% distinct( . )

times <- times  %>% arrange( PrecipGroup, period, plot, date ) %>% mutate ( )

head(times )

times$old_time <- times$date 

times 

head( times %>% group_by(PrecipGroup, period) %>% spread( plot , date ) ) 

times_list <- split(times, factor( paste0(times$PrecipGroup, times$period )) )

head( times_list[[4]] ) 


head( times_list[[4]] %>% ungroup(.) %>% spread(plot, date ) ) 


df <- expand.grid ( x = 1:10, y = 1:100)

df %>% spread( x, y )

head(times)
unique(times$id)
?spread



# Filter out bad readings ------------------------------------------------------------------------




# Correct switched ports -------------------------------------------------------------------------

stocks <- data.frame(
  time = as.Date('2009-01-01') + 0:9,
  X = rnorm(10, 0, 1),
  Y = rnorm(10, 0, 2),
  Z = rnorm(10, 0, 4)
)

stocksm <- stocks %>% gather(stock, price, -time)
head( stocksm ) 
stocksm %>% spread(stock, price)
stocksm %>% spread(time, price)
