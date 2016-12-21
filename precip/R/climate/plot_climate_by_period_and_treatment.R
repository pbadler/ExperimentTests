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
load('analysis/figure_scripts/my_plotting_theme.Rdata')
#------------------------------------------------------------------------ 
climate$Period[ climate$year <= 1925 ] <- 'not monitored'

climate <- climate %>% mutate( Period = ifelse( year < 2011, 'Historical', 'Modern'))
climate <- climate %>% filter( !(Treatment != 'Control' & year < 2011))
climate <- climate %>% unite( var, season, var) %>% spread(var, val  )


climate_long <- 
  climate %>% 
  filter( year < 2016) %>% 
  rename( `Fall Temperature` = fall_TAVG_avg, `Spring Temperature` = spring_TAVG_avg, `Summer Temperature` = summer_TAVG_avg, `Winter Temperature` = winter_TAVG_avg) %>%
  rename( `Mean Annual Temperature` = MAT) %>% 
  rename( `Winter PPT` = winter_PRCP_ttl, `Spring PPT` = spring_PRCP_ttl, `Summer PPT` = summer_PRCP_ttl, `Fall PPT` = fall_PRCP_ttl) %>% 
  mutate( `Total Annual PPT` = `Winter PPT` + `Spring PPT` + `Summer PPT` + `Fall PPT`) %>%
  select( year, Period, Treatment,  `Mean Annual Temperature`, `Total Annual PPT`, `Winter PPT`, `Spring PPT`, `Summer PPT`, `Fall PPT`, `Spring Temperature`, `Summer Temperature`, `Fall Temperature`, `Winter Temperature`) %>% 
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
  geom_ribbon( data = subset( x , Period == 'Historical'), aes( x = as.numeric(year), ymax = ymax, ymin = ymin), fill = 'gray', alpha = 0.01, color = 'gray')  + 
  geom_ribbon( data = subset( x , Period == 'Modern'), aes( x = as.numeric(year), ymax = ymax, ymin = ymin), fill = 'orange', alpha = 0.01, color = 'orange') + 
  geom_label(aes( x = as.numeric(label_x), y = label_y, label = period_label)) + 
  xlim( 1924, 2019) + 
  ylab(label =  parse( text = x$ylabel[1])) + 
    ggtitle( label = x$season[1]) + 
  scale_color_manual(values = c('black', 'red', 'blue')) + 
    my_theme
  
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
ts_plots$p[1]

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
      scale_fill_manual(values = c('gray', 'red')) +  
      guides(fill = guide_legend('Period', override.aes = list(size = 0, color = NA)), colour = 'none', shape = guide_legend('Treatment') )  + 
    ggtitle( label = x$season[1]) + 
    my_theme
    
}

rank_plots ( subset( rank_df , season == 'Total Annual PPT') ) 
rank_plots ( subset( rank_df, season == 'Mean Annual Temperature'))

rank_gg <- rank_df %>% group_by( season ) %>% do( p = rank_plots( . ))

pdf('figures/annual_climate_comparison_rank.pdf', height = 8 , width = 10)
print( rank_gg$p ) 
dev.off()

