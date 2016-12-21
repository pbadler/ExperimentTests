
rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(lme4)
library(zoo)

df <- readRDS('~/driversdata/data/idaho_modern/soil_moisture_data/data/processed_data/decagon_data_with_station_data.RDS')
load('analysis/figure_scripts/my_plotting_theme.Rdata')
# 
df <- df %>% mutate( v = ifelse(measure == 'VWC', v*100, v)) # convert to percent 
df$year <- as.numeric( strftime( as.Date( df$simple_date), '%Y' ) )
# summarize treatment differences:  -----------------------------------------------------------------------------------

plot_T_vals <- df %>% 
  filter( !is.na(v), bad_values == 0, measure == 'C', depth_label != '25 cm deep') %>%
  group_by(Treatment_label, season_label, tod, depth_label) %>% 
  summarise( avg = mean(v), stddev = sd(v), max = max(v), min = min(v), n = n(), ci = 1.96*(stddev/sqrt(n)), uci = avg + ci, lci = avg - ci ) 


air_T_diffs <- 
  df %>% 
  filter( measure == 'C', depth == 'air temperature', stat == 'raw', bad_values == 0) %>% 
  select( PrecipGroup, simple_date, season, tod, Treatment, v ) %>% 
  group_by( PrecipGroup, Treatment, season, simple_date, tod ) %>%
  summarise( avg_T = mean(v, na.rm = TRUE), n = n()) %>% 
  filter( !is.na(avg_T), n == 6)  


ggplot(subset( air_T_diffs, PrecipGroup == 3), aes( x = tod, y = avg_T, fill = Treatment )) + 
  geom_boxplot( ) + 
  facet_grid(  ~ season )+ 
  scale_fill_manual(values = my_colors[2:4])


ggplot(subset( air_T_diffs, PrecipGroup != 3), aes( x = tod, y = avg_T, fill = Treatment )) + 
  geom_boxplot( ) + 
  facet_grid(  ~ season )+ 
  scale_fill_manual(values = my_colors[2:4])


air_T_diffs <- air_T_diffs %>% spread(Treatment, avg_T) %>% mutate( DE = Drought - Control, IE = Irrigation - Control ) %>% dplyr::select( PrecipGroup, season, simple_date, tod, DE, IE ) %>% gather( Treatment, diff, DE, IE )

ggplot(air_T_diffs, aes( x = tod, y = diff, fill = Treatment ) ) + geom_boxplot() + facet_grid( PrecipGroup ~ season ) 
