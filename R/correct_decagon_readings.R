rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(zoo)

df <- readRDS(file = 'data/temp_data/decagon_data_corrected_dates.RDS')

df <- df %>% ungroup() 

# Filter out bad readings ------------------------------------------------------------------------

get_run_lengths <- function( x ) { 
  unlist ( lapply( rle(x)$lengths, function(x) rep(x, x) ) ) 
}

# My complicated steps for data classification:  ---------------------------------------------------------------
#   calculate rolling mean at 5 values 
#   calculate standard deviations from rolling mean 
#   calculate rolling mean of standard deviations 
#   classify values as high variability or low variabilty based on rolling mean of standard deviations 
#   classify sections as bad windows when 40/100 values are high variability 
#   calculate run lengths of good windows
#   only keep good windows of run lengths > 50 values in a row 

df<- df %>% 
  group_by(plot, port, measure) %>% 
  arrange(plot, port, measure, new_date ) %>% 
  mutate( frame_length = n() ) %>% 
  filter( frame_length > 8 ) %>% 
  mutate( rllm = rollapply( value, 8, mean, na.rm = TRUE, fill = NA, align = 'center')) %>% 
  mutate( dff = (value - rllm)^2) %>% 
  mutate( rllsd = rollapply( dff, 8, mean, na.rm = TRUE, fill = NA, align = 'center')) %>%
  mutate( avg_rllsd = rollapply( rllsd, 12, mean, na.rm = TRUE, fill = NA, align = 'center')) %>% 
  mutate( highv = ifelse( measure == 'VWC' & !is.na(avg_rllsd) & avg_rllsd > 0.001, 1, 0)) %>% 
  mutate( highv = ifelse( measure == 'C' & !is.na(avg_rllsd) & avg_rllsd > 400, 1, highv )) %>% 
  mutate( bad_window = rollapply( highv, 100 , sum, na.rm = TRUE, fill = NA)) %>% 
  mutate( bad_window = ifelse( !is.na(bad_window) & bad_window > 50, 1, 0)) %>%
  mutate( bad_window = ifelse( bad_window == 0 & measure == 'C' & (is.na(value) | value < -30 | value > 65 ) , 1, bad_window )) %>% 
  mutate( bad_window = ifelse( bad_window == 0 & measure == 'VWC' & (is.na(value) | value < -0.2 | value > 0.75), 1, bad_window )) %>% 
  mutate( bad_window = ifelse( bad_window == 0 | is.na(bad_window), 0, bad_window)) %>% 
  mutate( window_lengths = get_run_lengths( bad_window ) ) %>% 
  mutate( bad_window = ifelse( bad_window == 0 & window_lengths < 50, 1, bad_window)) %>% 
  gather( stat, v, value, rllm, rllsd, avg_rllsd) 

# manually remove 7_8_C port 4 last few months ---------------------------------------------------------------------------------------- 
df$bad_window <- as.numeric(df$bad_window)

df <- df %>% mutate(bad_window = ifelse( plot == '7_8_C' & port == 'Port 4' & measure == 'VWC' & new_date > strptime( '2015-07-01', format = '%Y-%m-%d'), 1, bad_window ) )

# --------------------------------------------------------------------------------------------------------------------------------------

df$stat <- factor(df$stat, label = c('avg rolling sd', 'rolling avg', 'rolling sd' , 'raw'))
df$plot <- as.character(df$plot)
df$bad_values <- factor(df$bad_window)

df <- df %>% 
  group_by(plot, port, period, measure ) %>% 
  #filter( good_date ==  1) %>% 
  mutate( has_vals = sum(stat == 'raw' & !is.na(v) ) > 0 ) %>%
  filter( has_vals)

p <- ggplot( df, aes ( x = new_date, y = v, color = stat))  + 
  geom_point() 

p2 <- ggplot( df , aes ( x = new_date, y = v, color = bad_values, group = depth)) + 
  geom_point()

pp <- df %>% 
  group_by(plot, port, period, measure) %>% 
  do( pp =  p %+%  . + ylab ( unique(.$measure)) + ggtitle( paste( c('plot ', 'period '), c(unique(.$plot), unique(.$period) ), collapse = '; ') ) ) 

pp2 <- df %>% 
  filter( stat == 'raw') %>%
  group_by( plot, port, period, measure) %>% 
  do(pp = p2 %+% . + ylab ( unique(.$measure)) + ggtitle( paste( c('plot ', 'period '), c(unique(.$plot), unique(.$period) ), collapse = '; ') ) ) 

pp_all <- df %>%
  filter( stat == 'raw' ) %>% 
  group_by( plot, port, measure) %>% 
  do(pp = p2 %+% . + ylab( unique( .$measure)) + ggtitle( paste( c( 'treatment:', unique(levels(.$Treatment)[.$Treatment]), '; plot:', unique(.$plot), ';', unique(.$port)), collapse = ' ') ) )

pdf( 'figures/check_bad_windows_continuous.pdf', height = 8, width = 10 )
print( pp_all$pp ) 
dev.off() 

#pdf( 'figures/check_bad_windows.pdf', height = 8, width = 10)
#print( pp2$pp ) 
#dev.off() 

# plot all readings 

plot_ts <- function( x ) { 
  
  ggplot ( data = x, aes( x = new_date, y = v , group = port, color = Treatment )) + 
    geom_point( alpha = 0.2, size = 0.2)   + 
    facet_wrap( ~ depth, ncol = 1 ) + 
    ylab( unique(x$measure)) + 
    ggtitle( paste( 'Plot Group:', unique(x$PrecipGroup)) )
  
}

plot_ts_diff <- function( x ) { 
  ggplot ( data = x, aes( x = new_date, y = v , color = Treatment )) + 
    geom_point( alpha = 0.2, size = 0.2)   + 
    facet_wrap( ~ depth, ncol = 1 ) + 
    ylab( paste( unique(x$measure), 'difference from control' )) + 
    ggtitle( paste( 'Plot Group:', unique(x$PrecipGroup)) )
}


df$depth <- factor(df$depth, labels = c('25 cm deep', '5 cm deep', 'air temperature'))   

p <- df %>% 
  filter( stat == 'raw', bad_values == 0 ) %>% 
  group_by(PrecipGroup, measure) %>% 
  do( p = plot_ts(x = . )  )


p_diff <- df %>% 
  filter( stat == 'raw', bad_values == 0 ) %>% 
  group_by( PrecipGroup, Treatment, measure, depth, new_date ) %>% 
  summarise( mean_v = mean( v ) ) %>% 
  select( PrecipGroup, Treatment, measure, depth, new_date, mean_v) %>% 
  spread( Treatment, mean_v) %>% 
  mutate ( Drought = Drought - Control, Irrigation = Irrigation - Control ) %>% 
  gather( Treatment, v, Control, Drought, Irrigation) %>%  
  group_by(PrecipGroup, measure) %>% 
  filter( Treatment != 'Control') %>%
  do( p = plot_ts_diff(x = . )  )


pdf('figures/corrected_decagon_time_series.pdf', width = 10, height = 7)
print( p$p )   
dev.off() 

pdf( 'figures/corrected_decagon_time_series_plot_differences.pdf', width = 10, height = 7) 
print ( p_diff$p ) 
dev.off() 


saveRDS(df, 'data/temp_data/decagon_data_corrected_values.RDS')


