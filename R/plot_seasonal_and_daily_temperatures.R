# check date helper functions 

quick_plot_check <- function( x ) { 
  
  plot_check <- x %>% ungroup ( ) %>% filter(depth == 'air')
  par(mfrow = c(1,2))
  plot( data= plot_check, value ~ new_date, ylab = 'air', main = 'temp by date')
  plot( data= plot_check, value ~ new_hour, ylab = 'air', main = 'temp by hour')
  par(mfrow = c(1,1))
} 


out %>% group_by ( f ) %>% do ( quick_plot_check ( . ))

quick_plot_check(x = out %>% filter( f == 'data/soil_moist_data/2016_Spring/EL5739 10May16-1316.txt'))

