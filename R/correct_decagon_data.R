rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)

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

# Filter out bad readings ------------------------------------------------------------------------

port_stats <- df %>% 
  group_by( plot, port, period, measure ) %>% 
  summarise( mean = mean(value, na.rm = TRUE), 
             sd  = sd( value, na.rm = TRUE ), 
             mx = max(value, na.rm = TRUE), 
             mn = min(value, na.rm = TRUE), 
             range = mx - mn )

vwc_stats <- port_stats %>%
  ungroup( ) %>% 
  filter( measure == 'VWC') %>%
  arrange(  desc( range  ), desc(sd) )  

# loop through each period and visually inspect the VWC data 

# p <- ggplot( data = subset(df, plot == 1 & port == 'Port 2' & period == 1), aes( x = date, y = value)) + 
#   geom_point() + 
#   ylim ( -1, 1) 
# 
# pdf( 'figures/VWC_series_by_port_and_period.pdf', height = 8, width = 10)
# for ( i in 1:nrow( subset( vwc_stats, !is.na(sd)))) { 
#   print( p %+% subset( df , measure == 'VWC' & plot == vwc_stats$plot[i] & port == vwc_stats$port[i] & period == vwc_stats$period[i])  )
#   }
# dev.off()

# suggest using a lagged difference to filter out points that jump around too much
library(zoo)

get_run_lengths <- function( x ) { unlist ( lapply( rle(x)$lengths, function(x) rep(x, x) ) ) }  

df<- df %>% 
  group_by(plot, period, port, measure) %>% 
  arrange(plot, port, measure, date ) %>% 
  mutate( rllm = rollapply( value, 5, mean, na.rm = TRUE, fill = NA, align = 'center')) %>% 
  mutate( dff = (value - rllm)^2) %>% 
  mutate( rllsd = rollapply( dff, 5, mean, na.rm = TRUE, fill = NA, align = 'center')) %>%
  mutate( avg_rllsd = rollapply( rllsd, 12, mean, na.rm = TRUE, fill = NA, align = 'center')) %>% 
  mutate( highv = ifelse( measure == 'VWC' & !is.na(avg_rllsd) & avg_rllsd > 0.0001, 1, 0)) %>% 
  mutate( highv = ifelse( measure == 'C' & !is.na(avg_rllsd) & avg_rllsd > 20, 1, highv )) %>% 
  mutate( bad_window = rollapply( highv, 100 , sum, na.rm = TRUE, fill = NA)) %>% 
  mutate( bad_window = ifelse( !is.na(bad_window) & bad_window > 40, 1, 0)) %>%
  mutate( bad_window = ifelse( bad_window == 0 & measure == 'C' & (is.na(value) | value < -30 | value > 55 ) , 1, bad_window )) %>% 
  mutate( bad_window = ifelse( bad_window == 0 & measure == 'VWC' & (is.na(value) | value < -0.2 | value > 0.75), 1, bad_window )) %>% 
  mutate( window_lengths = get_run_lengths( bad_window ) ) %>% 
  mutate( bad_window = ifelse( bad_window == 0 & window_lengths < 50, 1, bad_window)) %>% 
  gather( stat, v, value, rllm, rllsd, avg_rllsd) 

df$stat <- factor(df$stat, label = c('avg rolling sd', 'rolling avg', 'rolling sd' , 'raw'))
df$plot <- as.character(df$plot)
df$bad_values <- factor(df$bad_window)

test_case <- df %>% 
  group_by(plot, port, period, measure ) %>% 
  filter( good_date ==  1) %>% 
  mutate( has_vals = sum(stat == 'raw' & !is.na(v) ) > 0 ) %>%
  filter( has_vals)
  
p <- ggplot( test_case, aes ( x = date, y = v, color = stat))  + 
  geom_point() 

p2 <- ggplot( test_case , aes ( x = date, y = v, color = bad_values, group = depth)) + 
  geom_point()

pp <- test_case %>% 
  group_by(plot, port, period, measure) %>% 
  do( pp =  p %+%  . + ylab ( unique(.$measure)) + ggtitle( paste( c('plot ', 'period '), c(unique(.$plot), unique(.$period) ), collapse = '; ') ) ) 

pp2 <- test_case %>% 
  filter( stat == 'raw') %>%
  group_by( plot, port, period, measure) %>% 
  do(pp = p2 %+% . + ylab ( unique(.$measure)) + ggtitle( paste( c('plot ', 'period '), c(unique(.$plot), unique(.$period) ), collapse = '; ') ) ) 

pp_all <- test_case %>%
  filter( stat == 'raw' ) %>% 
  group_by( plot, measure) %>% 
  do(pp = p2 %+% . + ylab( unique( .$measure)) + ggtitle( paste( c( 'plot:', unique(.$plot), ';', unique(.$port)), collapse = ' ') ) )

pdf( 'figures/check_bad_windows_continuous.pdf', height = 8, width = 10 )
pp_all$pp
dev.off() 

pdf( 'figures/check_bad_windows.pdf', height = 8, width = 10)
pp2$pp
dev.off() 



#



# Correct switched ports -------------------------------------------------------------------------

