clim <- readRDS('data/temp_data/all_clim_covs.RDS')


head( clim )

clim <- subset( clim , year > 2010 )

clim <- clim %>% dplyr::select(Treatment, Period, year, quarter,  contains("sp.1"))

clim <- 
  clim %>% 
  gather( var, val , P.f.w.sp.1:VWC.sp.1)  %>% 
  filter( var != 'T.sp.1') 

my_colors <- c('#66c2a5','#fc8d62','#8da0cb')

ggplot( clim %>% spread( var, val ), aes( x = P.f.w.sp.1, y = VWC.sp.1, color = Treatment )) + 
  geom_point() + 
  xlab('Total fall winter-spring precip mm ') + 
  ylab('Modeled spring shallow soil moisture (%) ')  + 
  scale_color_manual(values = my_colors) + 
  theme_bw()
  

