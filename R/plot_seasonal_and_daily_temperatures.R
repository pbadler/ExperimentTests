rm(list = ls())
library( dplyr)
library(ggplot2)

temps <- readRDS('data/temp_data/decagon_data_with_station_data.RDS')

df  <- temps %>% 
  filter( bad_values == "0", depth == 'air temperature', !is.na(v) , stat == 'raw') 

test <- df %>% ungroup() %>% filter( f == f[1] )

daily_t <- ggplot(test, aes( x  = hour , y = v , group = simple_date ) ) +  
  geom_line(stat = 'smooth', method = 'loess', alpha = 0.1) + 
  geom_vline(aes(xintercept = 16)) + 
  ylim ( -50, 60 ) + 
  facet_wrap( ~ month )

plots <- df %>% ungroup() %>% 
  group_by (plot , year ) %>% 
  do ( gg = daily_t %+% . + ggtitle ( paste( 'plot = ', unique( .$plot) , 'year = ', unique( .$year) ))) 

pdf ( 'figures/check_daily_temperatures.pdf', width = 7 , height = 10)
print( plots$gg ) 
dev.off()

#

daily_Tmeans <- df %>% 
  group_by ( plot , simple_date) %>% 
  mutate( tmean = mean(v), n = n()) %>% 
  filter( n == 12 ) %>% 
  group_by ( simple_date ) %>%
  mutate( tmean_all = mean(v ), 
          doy = as.numeric(strftime(new_date, '%j'))) 

daily_T_plot <- ggplot(daily_Tmeans, aes( x = doy, y = tmean, color = plot  )) + 
  geom_point() + 
  geom_line() + facet_wrap(~ year, ncol = 1 ) 

p <- daily_Tmeans %>% group_by(PrecipGroup, year ) %>% do( plot = daily_T_plot %+% . + ggtitle(paste( .$PrecipGroup, .$year)))

p$plot


# 

test <- 
  temps %>% 
  ungroup() %>% 
  filter(  stat == 'raw', plot == 15, depth == 'air temperature', simple_date > strftime( '2016-02-01', '%Y-%m-%d', tz = 'MST'), simple_date < strftime( '2016-03-01', '%Y-%m-%d', tz = 'MST')) %>% 
  select( new_date, simple_date, date.x,  reading, plot,  v, tod ) %>%
  arrange( new_date ) 


ggplot( test, aes( x = new_date, y = v, color = tod , group = simple_date )) + geom_point() + geom_line() 

#


annual_t <- ggplot ( test, aes( x = new_date, y = v , group = plot )) + geom_point(alpha = 0.1 ) + ylim ( -20, 60 ) 

plots <- df %>% ungroup () %>% group_by( plot) %>% do ( gg = annual_t %+% . + ggtitle (paste( 'plot = ' , unique(.$plot))))

plots$gg 


df_soil  <- temps %>% 
  filter( bad_values == "0", depth == '5 cm deep', measure == 'VWC' , !is.na(v) , stat == 'raw') %>% 
  mutate( simple_date = as.Date ( new_date ), year = strftime( new_date, '%Y'), month = strftime( new_date, '%m'))

annual_soil <- ggplot ( test, aes(x = new_date, y = v , group = plot )) + geom_point( alpha = 0.01 ) + ylim ( -0.2, 0.6)

plots <- df_soil %>% ungroup () %>% group_by( plot) %>% do ( gg = annual_soil %+% . + ggtitle (paste( 'plot = ' , unique(.$plot))))
plots$gg

