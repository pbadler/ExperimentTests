# save plotting theme to be used on all plots 


my_theme <- 
  theme_bw () + 
  theme ( panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank())


# black historical 
# olive control-modern 
# orange drought 
# blue  irrigation 

my_colors <- c('black', '#1b9e77', '#d95f02', '#7570b3') 


pdf_settings <- c('height' = 5, width = 5 ) 


save(my_colors, my_theme, file = 'analysis/figure_scripts/my_plotting_theme.Rdata')
