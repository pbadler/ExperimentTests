rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(zoo)

df <- readRDS(file = 'data/temp_data/decagon_data.RDS')

# correct bad dates  ------------------------------------------------------------------------------ 

find_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

fill_in_hours_skipped <- function( x ) { 
  
  hs = 0 
  
  for( i in 1:nrow(x)) {
    
    if (is.na( x$change[i] )) {
      
      x$hours_skipped[i] <- hs
      
    }else if(x$change[i] == 1 ){
      
      print(paste('old hs', hs ))
      
      hs <- x$hours_skipped[i] <- x$hours_skipped[i] + hs
      
      print(paste('new hs', hs))
      
    }else if(x$change[i] == 0 ){
      
      hs <- x$hours_skipped[i] <- 0 } 
  }
  
  return( x )
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

df <- df %>% filter( f != 'data/soil_moist_data/2015_Spring/EM20068_2015-04-30-0957.txt')

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
  filter( f != 'data/soil_moist_data/2015_Spring/EM20068_2015-04-30-0957.txt') %>% 
  arrange( date, f  ) 

# determined for each jump whether it should be corrected or remain in place 
# change = 1  indicates jumps that should be changed 
change <- 
  c( 1, 1, 1, 1, 1, 1, 1, 
   0, 0, 0, 0, 0, 0, 0, 0, 
   1, 1, 
   0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0, 
   1, 
   0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0 ,0, 0, 
   1, 1)

check$change  <- change 

# need to mark the initial row of this file as a change manually ------------------------
date1 <- df %>%  filter( f == 'data/soil_moist_data/2014_Fall/11_15Sep14-1708.txt', row_number() == 1 ) %>% ungroup() %>% select(date) %>% distinct( date)
old_date <- date1$date
new_date <- strptime( "2014-04-13 12:00:00", format = "%Y-%m-%d %H:%M:%S", tz = 'MST') 

add_time_diff <- as.numeric( diff.POSIXt(c(new_date, old_date)), tz = 'MST' , units = 'secs')

add_row <- df %>% filter( f == 'data/soil_moist_data/2014_Fall/11_15Sep14-1708.txt') %>% ungroup() %>% filter( row_number() == 1 ) %>% select(f, date, reading)  
add_row <- add_row %>% mutate( time_diff = add_time_diff,  hours_skipped = time_diff/(60*60) , reading_diff = 1, jump = 1, change = 1)
check <- rbind( check, add_row )

#-----------------------------------------------------------------------------------------
write.csv(check, 'data/temp_data/check_dates.csv', row.names = FALSE) # write list of changes 
# ------------------------------------------------------------------------------------------- 

df <- left_join(df, check , by =c( 'f', 'date', 'reading' )) # join changes to main df 

df <- df %>% 
  ungroup () %>% 
  group_by(f, plot, port, measure ) %>% 
  arrange( reading ) %>% 
  mutate( hours_skipped = ifelse( row_number() == 1 & is.na(change), 0, hours_skipped ))

out <- df %>%  do ( fill_in_hours_skipped(. ) ) # apply fill in hours function to all measurement groups 

# actually make the date changes here ----------------------------------------------------------------------------------

out <- out %>% 
  mutate( new_date = as.POSIXct(date - 60*60*hours_skipped, origin = '1970-01-01 00:00:00', tz = 'MST'))

# ----------------------------------------------------------------------------------------------------------------------
out <- out %>% 
  mutate ( good_date = ifelse ( date >= lag_mod_date - 60*60*12 & date <= modified_date + 60*60*12 , 1, 0))

# only take first record of duplicates  ------------------------------------
out <- out %>% filter( !is.na(value) ) %>% group_by (plot, port, depth, measure, new_date, value ) %>% arrange( f ) %>% filter( row_number() == 1 ) 

# check earliest and latest dates -----------------------------------------------------------------
out %>% ungroup( ) %>% summarise ( max( new_date ), min( new_date ), which.min(new_date ), which.max(new_date ))

# ---------------------------------------------------------------------------- 

out <- out %>% 
  ungroup() %>%
  mutate( simple_date = as.Date(new_date, tz = 'MST'), 
          hour = strftime( new_date, '%H', tz = 'MST'), 
          year = strftime( new_date, '%Y', tz = 'MST'), 
          month = strftime( new_date, '%m', tz = 'MST'))

saveRDS( out , 'data/temp_data/decagon_data_corrected_dates.RDS')

# ----------------------------------------------------------------------------