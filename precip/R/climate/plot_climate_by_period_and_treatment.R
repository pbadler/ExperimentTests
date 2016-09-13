###########################################################################
#
# Compare climate of historical Period to long-term average 
# and compare climate of contemporary Period to long-term average
#
###########################################################################

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

climate <- readRDS('data/temp_data/seasonal_climate.RDS')

#------------------------------------------------------------------------ 

my_colors <- c('#fc8d62', 'gray', '#8da0cb')
climate$Period[ climate$year == 1925 ] <- 'not monitored'

climate_long <- 
  climate %>% 
  filter( year < 2016) %>% 
  rename( `Fall Temperature` = fall_TMEAN_l0, `Spring Temperature` = spring_TMEAN_l0, `Summer Temperature` = summer_TMEAN_l0, `Winter Temperature` = winter_TMEAN_l0) %>%
  rename( `Mean Annual Temperature` = MAT) %>% 
  rename( `Cool Season PPT` = cool_PRCP_l0, `Warm Season PPT` = warm_PRCP_l0) %>% 
  mutate( `Total Annual PPT` = `Cool Season PPT` + `Warm Season PPT`) %>%
  select( year, Period, Treatment,  `Mean Annual Temperature`, `Total Annual PPT`, `Cool Season PPT`, `Warm Season PPT`, `Spring Temperature`, `Summer Temperature`, `Fall Temperature`, `Winter Temperature`) %>% 
  gather( season, value , `Mean Annual Temperature`:`Winter Temperature`) %>% 
  mutate( ylabel = ifelse( str_detect(season, pattern = 'PPT'), 'Total~precipitation~(mm)', 'Average ~ air~temperature~degree*C')) 


annual_plots <- 
  climate_long %>% 
  #filter(Treatment == 'Control') %>%
  group_by(Period) %>% 
  mutate( label_x = median(year)) %>%
  mutate(period_label = paste0(Period, '\n', min(year), '-', max(year))) %>%
  mutate(period_label = ifelse(Period == 'not monitored', 'not monitored', period_label)) %>% 
  ungroup() %>% 
  group_by(season) %>% 
  mutate( ymax = max(value) + 0.05*(max(value) - min(value)), 
          label_y = min(value) - 0.1*(max(value)- min(value)), 
          ymin = min(value) - 0.2*(max(value)- min(value))) 

base_ts <- function( x ) { 
  g <- ggplot( x, aes( x = year, y = value, group = 1)) +
  geom_ribbon( data = subset( x , Period == 'Historical'), aes( x = as.numeric(year), ymax = ymax, ymin = ymin), fill = my_colors[3], alpha = 0.3, color = my_colors[3])  + 
  geom_ribbon( data = subset( x , Period == 'Modern'), aes( x = as.numeric(year), ymax = ymax, ymin = ymin), fill = my_colors[1], alpha = 0.3, color = my_colors[1]) + 
  geom_label(aes( x = as.numeric(label_x), y = label_y, label = period_label)) + 
  xlim( 1924, 2019) + 
  ylab(label =  parse( text = x$ylabel[1])) + 
    ggtitle( label = x$season[1]) + 
  scale_color_manual(values = c('black', 'red', 'blue')) 
  
  if(str_detect(x$season[1], "PPT")){ 
    g +   
      geom_point(aes(shape = Treatment, color = Treatment)) + 
      geom_line(aes( group = Treatment, linetype = Treatment, color = Treatment))  
  }else{ 
    g + 
      geom_point() + 
      geom_line() 
  }
}

ts_plots <- annual_plots %>% group_by( season ) %>% do( p = base_ts(.) )

pdf('figures/annual_climate_comparison_timeseries.pdf', height = 8 , width = 10)
print( ts_plots$p ) 
dev.off()

#------------------------------------------------------------------------------- 

rank_df <- 
  annual_plots %>% 
  group_by(season) %>% 
  filter( !( !str_detect(pattern = 'PPT', season) & Treatment != 'Control') ) %>% 
  mutate( rank_order = rank(`value`, ties.method = 'first')) %>% 
  mutate( xlab = ifelse( str_detect(pattern = "PPT", season ), 'Year rank from driest to wettest', 'Year rank from coldest to hotest'))

rank_plots <- function( x ) { 

  ggplot( data = x, aes( x = rank_order, y = value, fill = Period)) + 
      geom_bar(stat = 'identity') + 
      geom_point( data = subset(x , Treatment != 'Control'), aes( y = 1.05*value, shape = Treatment, color = NULL), size = 6) + 
      ylab(label =  parse( text = x$ylabel[1])) +
      xlab(label =  x$xlab[1]) + 
      scale_shape_manual(values = c('-', '+')) + 
      scale_fill_manual(values = my_colors) +  
      guides(fill = guide_legend('Period', override.aes = list(size = 0, color = NA)), colour = 'none', shape = guide_legend('Treatment') )  + 
    ggtitle( label = x$season[1]) 
    
}

rank_plots ( subset( rank_df , season == 'Total Annual PPT') ) 
rank_plots ( subset( rank_df, season == 'Mean Annual Temperature'))

rank_gg <- rank_df %>% group_by( season ) %>% do( p = rank_plots( . ))

pdf('figures/annual_climate_comparison_rank.pdf', height = 8 , width = 10)
print( rank_gg$p ) 
dev.off()

