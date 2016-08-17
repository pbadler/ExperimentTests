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

df <- df %>% mutate( hour = strftime( date , '%H'))

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

reading_list <- df %>% ungroup () %>% select( f, plot, id , period, date, reading ) %>% mutate( f = factor(f)) %>% distinct()

table( reading_list$f, reading_list$period ) # one file per period 

jumps <- reading_list %>% 
  group_by(f ) %>% 
  arrange( f , reading ) %>% 
  mutate( time_numeric = as.numeric(date )) %>% 
  mutate ( time_diff = c(NA, diff(time_numeric, 1 ))) %>%
  mutate( hours_skipped = time_diff/3600 - 2 ) %>% 
  mutate( reading_diff = c(NA, diff(reading, 1))) %>% 
  ungroup() %>% 
  mutate( jump = ifelse( reading_diff == 1 & (hours_skipped != 0 ), 1, 0 )) %>% 
  mutate( lead_jump = lead( jump, 1 )) 


jumps %>% group_by ( f ) %>% summarise( n_jumps =  sum(jump, na.rm = T)) %>% filter ( n_jumps > 0  )
  
check <- 
  jumps %>% 
  filter ( f != 'data/soil_moist_data/2014_Fall/15_15Sep14-1802.txt' ) %>% 
  select( f, date, reading , time_diff, hours_skipped, reading_diff, jump ) %>% 
  filter( jump > 0 , (hours_skipped != 0 ) & reading_diff == 1 ) %>% 
  filter( ! ( hours_skipped == 2 & f == 'data/soil_moist_data/2015_Fall/EM6845 4Nov15-1123.txt')) %>% 
  filter( ! (hours_skipped == 2 & f == 'data/soil_moist_data/2016_Spring/EM6845 10May16-1458.txt')) %>% 
  filter( f != 'data/soil_moist_data/2015_Fall/EL5739 4Nov15-1838.txt') %>% 
  filter( f != 'data/soil_moist_data/2015_Fall/EL5742 4Nov15-1820.txt') %>% 
  filter( !( abs(hours_skipped) < 10000 & f == 'data/soil_moist_data/2015_Fall/EL5743 4Nov15-1828.txt')) %>% 
  filter( f != 'data/soil_moist_data/2013_Spring/EM20070.txt') %>% 
  filter( f != 'data/soil_moist_data/2013_Spring/EM20085.txt') %>% 
  arrange( date, f  ) 

# determined for each jump whether it should be corrected or remain in place 
# change = 1  indicates jumps that should be changed 
change <- 
  c( 1, 1, 1, 1, 1, 1, 1, 
   0, 0, 0, 0, 0, 0, 0, 0, 
   1, 1, 
   0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0, 
   1, 1, 
   0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0, 
   1, 1)

check$change  <- change 

write.csv(check, '~/Downloads/check_dates.csv', row.names = FALSE)

quick_plot_check <- function( x = df, file_name = 'data/soil_moist_data/2015_Fall/EL5739 4Nov15-1838.txt' ) { 
  
  plot_check <- x %>% ungroup ( ) %>% filter( f == file_name , depth == 'air')
  par(mfrow = c(1,2))
  plot( data= plot_check, value ~ date, ylab = 'air', main = 'temp by date')
  plot( data= plot_check, value ~ hour, ylab = 'air', main = 'temp by hour')
  par(mfrow = c(1,1))
} 

quick_plot_check(file_name = 'data/soil_moist_data/2015_Fall/EL5739 4Nov15-1838.txt')

#df <- left_join( df, check[ , c('f', 'date', 'reading', 'change', 'time_diff', 'hours_skipped')], by = c('f', 'date', 'reading'))

df <- df %>% arrange( f, reading )

df$time_adj <- 0 



check %>% arrange ( f, reading )

for( i in 1:nrow( check)  ){ 
  
  file = check$f[i]  
  rd = check$reading[i]
  adj = check$hours_skipped[i]
  df <- df %>% ungroup() %>% mutate( time_adj = ifelse (as.character( f )== as.character(file ) & reading >= rd , adj, time_adj ))
  
  df <- df %>% ungroup( ) %>% mutate( time_adj = ifelse ( as.character( f)  == as.character( file) & reading >= rd , adj, time_adj  ))
} 

df %>% mutate( adjusted )


df %>% 
  group_by( f ) %>% 
  arrange( f, desc( reading ) ) %>%
  
  mutate( change = ifelse ( row_number() == 1, 0, change ) ) %>%
  mutate( change = cumsum(change)) %>%  
  mutate( change = ifelse ( is.na(change), lag(change, 1), change ) ) %>% 
  dplyr::select( change )




df %>% arrange( f, desc(reading ) ) 

df <- df %>% 
  mutate ( good_date = ifelse ( date >= lag_mod_date - 60*60*12 & date <= modified_date + 60*60*12 , 1, 0))

saveRDS( df , 'data/temp_data/decagon_data_corrected_dates.RDS')
