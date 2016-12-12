# save plotting theme to be used on all plots 


my_theme <- theme_bw () + theme ( panel.grid.major = element_blank(), 
                                  panel.grid.minor = element_blank())

my_colors <- c('#1b9e77', '#d95f02', '#7570b3')

save(my_colors, my_theme, file = 'analysis/figure_scripts/my_plotting_theme.Rdata')
