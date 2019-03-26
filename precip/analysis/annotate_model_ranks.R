rm(list = ls())
library(tidyverse)

files <- dir('output', '.*_ranks3.csv$', full.names = T)
files <- c(files, dir('output', '.*_survival_model_ranks2.csv$',full.names = T))

ranks <- lapply( files, read_csv)

ranks <- do.call(rbind, ranks) 

ranks <- 
  ranks %>% 
  group_by(spp, vr ) %>% 
  gather(parameter, est, starts_with('beta')) 

labels <- data.frame( parameter = unique( ranks$parameter ), 
            description = c('intercept', 'size', 'small_plants', 'intra_comp', 'inter_comp', 'Temp', 'Moist', 'Temp_x_Moist')) 


ranks <- 
  ranks %>% 
  left_join(labels) %>% 
  select(-parameter) %>% 
  spread( description, est ) %>% 
  arrange( spp, vr, desc( out_of_sample_lppd) ) %>% 
  ungroup %>% 
  mutate( climate_window = str_extract(climate_effects, '\\d+')) %>% 
  mutate( climate_window = ifelse( is.na(climate_window), 'NULL_MOD', climate_window)) 
  

window_descriptions <- 
  ranks %>%
  distinct(climate_window, Temperature) %>% 
  mutate( seasons = str_remove_all(Temperature, 'C\\.T\\.')) %>% 
  arrange( climate_window)  %>% 
  select( climate_window, seasons)
  

ranks <- 
  ranks %>% 
  rename( 'oos_lppd' = out_of_sample_lppd, 
          'oos_mse' = out_of_sample_mse) %>% 
  group_by(spp, vr) %>% 
  mutate( model_rank = row_number()) %>% 
  select( vr, spp, model_rank, climate_window, oos_lppd, oos_mse, 
          intercept, size, small_plants, intra_comp, inter_comp, Temp, Moist, Temp_x_Moist) 
  

write_csv(path = '~/Dropbox/projects/ExperimentTests/precip/output/stan_model_ranks.csv', ranks)
#write_csv(path = '~/Desktop/climate_window_descriptions.csv', window_descriptions)
