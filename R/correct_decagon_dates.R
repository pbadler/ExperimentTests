rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(zoo)

df <- readRDS(file = 'data/temp_data/decagon_data.RDS')

# Filter out redundant and bad dates ------------------------------------------------------------------------------ 

find_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

lag_period <- df %>% 
  ungroup (.) %>% 
  select ( PrecipGroup, plot, period, modified_date ) %>%
  distinct() %>%
  mutate (period = period + 1 ) %>% 
  rename( lag_mod_date = modified_date)

df <- df %>%   
  arrange( PrecipGroup, plot, period ) %>% 
  group_by(PrecipGroup, plot, period ) %>%
  left_join( lag_period, by = c('PrecipGroup', 'plot', 'period'))

df$date_started <- as.character ( levels( df$date_started )[df$date_started] )  

df <- df %>% 
  mutate( lag_mod_date = ifelse(is.na(lag_mod_date), as.POSIXct(date_started), lag_mod_date)) %>% 
  mutate( lag_mod_date = as.POSIXct(lag_mod_date, origin = '1970-01-01', tz = 'MST'))

df <- df %>% 
  mutate ( good_date = ifelse ( date >= lag_mod_date - 60*60*12 & date <= modified_date + 60*60*12 , 1, 0))

saveRDS( df , 'data/temp_data/decagon_data_corrected_dates.RDS')
