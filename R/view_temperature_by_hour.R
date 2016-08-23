rm(list = ls () )
library( dplyr ) 
library( tidyr )

temps <- readRDS('data/temp_data/decagon_data_corrected_values.RDS')

temps <- temps %>% mutate( hour = as.numeric(strftime( date, '%H') ), doy = as.numeric(strftime(date, '%j')), month = as.numeric( strftime( date, '%m')) )

good_temps <- 
  temps %>% 
  filter(good_date == 1, bad_values == 0, measure == 'C' , depth == 'air temperature' ) %>% 
  group_by( plot, port , period, doy) %>% 
  mutate( max_v = max(v), min_v = min(v), sub_v = max_v - v ) %>% 
  ungroup() 

nrow(good_temps)
class(good_temps$date)
class(good_temps$date_started)



hist_temps <- good_temps %>%
  filter( doy > 60 ) %>% 
  filter( v == max_v )

hist_time <- ggplot(hist_temps, aes( x = hour) ) + facet_grid(PrecipGroup +  plot ~ . ) + geom_histogram() 

gg_out <- hist_temps %>% group_by( period ) %>% do( gg = hist_time %+% . + ggtitle(unique(.$period)) )

gg_out$gg

time_plot <- ggplot ( good_temps, aes( x = hour, y = v, group = doy)) + 
  geom_point( alpha = 0.1 ) + 
  #geom_line(stat = 'smooth', method = 'loess', alpha = 0.1) + 
  facet_grid(Treatment ~ PrecipGroup) + 
  ylim( c(-20, 40))

hour_plots <- good_temps %>% filter( month %in% c(8)) %>% group_by(period) %>% do ( gg = time_plot %+% . + ggtitle(unique(.$period)))
hour_plots$gg

good_temps_test  <- good_temps %>% mutate( hour  =  ifelse(PrecipGroup == 3 & period == 4, hour - 6,  hour ), hour = ifelse(hour < 0 , hour + 24, hour ))

hour_plots <- good_temps_test %>% filter( month %in% c(8)) %>% group_by(period) %>% do ( gg = time_plot %+% . + ggtitle(unique(.$period)))
hour_plots$gg
