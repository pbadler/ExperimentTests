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





# Filter out bad readings ------------------------------------------------------------------------




# Correct switched ports -------------------------------------------------------------------------