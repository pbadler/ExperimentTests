rm(list = ls())
library(tidyverse)

#lmer_mods <- read_csv('~/Dropbox/projects/ExperimentTests/precip/output/lmer_growth_model_ranks.csv')

stan_mods <- read_csv('~/Dropbox/projects/ExperimentTests/precip/output/stan_model_ranks.csv')

top_mods <- 
  stan_mods %>% 
  group_by( vr, spp) %>% 
  filter( row_number() == 1 ) %>% 
  select( vr, spp, climate_window, Temp:Temp_x_Moist) 

gg1 <- 
  top_mods %>% 
  gather( par, est, Temp:Temp_x_Moist) %>% 
  mutate( col = ifelse(est < 0, 'red', 'blue')) %>%  
  ggplot(aes( x = par, y = est, color = col)) + 
  geom_point() + 
  geom_hline(aes( yintercept = 0), linetype = 2) + 
  facet_grid(vr~spp, scales = 'free_y') + 
  scale_color_manual(values = c('blue', 'red'))

plot_2d_grid <- list()
obs_years <- list()

for( i in 1:nrow(top_mods)){ 
  vr <- top_mods$vr[i]
  spp <- top_mods$spp[i]
  
  fn <- paste0('data/temp_data/', spp, '_growth_survival_dataframe.RDS')

  temp <- readRDS(fn)
  
  top_mod <- top_mods$climate_window[top_mods$spp == spp & top_mods$vr == vr ]
  
  if(top_mod == "NULL_MOD"){ 
    top_mod <- 6
  }
  
  top_covs <- paste0(c('C.VWC.', 'C.T.'), top_mod)
  
  all_years <- 
    temp %>%
    filter( Period == 'Historical') %>% 
    select(year, top_covs) %>% 
    distinct() 
  
  names( all_years ) <- c('year', 'moist', 'therm')
  
  vars <- expand.grid( temp_x = seq(-2, 3.5, by = 0.1), 
                       moist_y = seq(-2, 2, by = 0.1))
  
  plot_2d_grid[[i]] <- 
    data.frame(vars, top_mods[i, ]) %>% 
    mutate( response = temp_x*Temp + moist_y*Moist + temp_x*moist_y*Temp_x_Moist ) %>% 
    select( spp, vr, temp_x, moist_y, response )
  
  all_years$spp <- spp 
  all_years$vr <- vr 
  obs_years[[i]] <- all_years 
  if( all(is.na(plot_2d_grid[[i]]$response))){ 
    obs_years[[i]]$therm <- NA
    obs_years[[i]]$moist <- NA
  }
  
}

plot_2d_grid <- do.call(rbind, plot_2d_grid)
obs_years <- do.call(rbind, obs_years)

gg2 <- 
  plot_2d_grid %>% 
  group_by(vr) %>% 
  mutate( scaled_response = scale(response)) %>% 
  ggplot( aes( x = moist_y, y = temp_x, z = scaled_response, fill = scaled_response )) + 
  geom_tile() + 
  geom_contour(bins = 4) + 
  geom_point(data = obs_years, aes( x = moist, y = therm, z = NULL, fill = NULL ), alpha = 0.5) + 
  facet_grid(vr~spp) + 
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', name = 'response predicted') + 
  xlab( 'Annual Soil Moisture\nDry      --->       Wet') + 
  ylab( 'Annual Temperature\nCold      --->      Hot') + 
  ggtitle('Predicted effects of temperature and moisture in top models')

gg2 + coord_equal()

ggsave(gg2 + coord_equal(), filename = 'figures/top_model_responses.png', width = 8, height = 6)

ggsave( 
  gg1 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + 
  xlab('parameter') + 
  ylab('estimate') + 
  ggtitle( 'climate effects in top models'), 
  filename = 'figures/top_model_effects.png', width = 8, height = 6)

