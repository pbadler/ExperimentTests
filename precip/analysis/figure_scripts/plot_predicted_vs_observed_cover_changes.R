# 
library(dplyr)
library(tidyr)
library(stringr)

rm(list = ls() ) 
spp = 'ARTR'

# get last year of cover which is not in the survival dataframe ------------------------------------- #   
oldCover <- read.csv(paste0( '~/driversdata/data/idaho/speciesData/', spp, '/quadratCover.csv'))
newCover <- read.csv(paste0( '~/driversdata/data/idaho_modern/speciesData/', spp, '/quadratCover.csv'))
oldCover$Period <- "Historical"
newCover$Period <- "Modern"
last_cover  <- rbind(oldCover, newCover)
quad <- read.csv('~/driversdata/data/idaho_modern/quad_info.csv')

last_cover <-
  last_cover %>%
  left_join(quad) %>%
  mutate( quad = as.numeric(str_extract(pattern = '[0-9]+$', quad))) %>% 
  filter( Treatment %in% c('Control', 'Drought', 'Irrigation')) %>%
  mutate( year = ifelse(year < 100, year + 1900, year)) %>% 
  mutate( observed = totCover)

# ------------------------------------------------------------------------------------------------------ #
df <- read.csv(paste0('data/temp_data/', spp, '_growth_and_survival_cleaned_dataframe.csv'))

df$true_cov <- exp(df$logarea.t0) # use year 0 cover 

trueCov <- 
  df %>% 
  group_by( Period, year, Treatment, quad ) %>% 
  summarise( observed = sum( true_cov ) )

obs_cover <- merge( last_cover, trueCov, by = c('Period', 'year', 'Treatment', 'quad'), all.x = T)

obs_cover <- obs_cover %>% mutate( observed = ifelse( is.na(observed.y), observed.x, observed.y)) # fill in data from end years in dataframe


# get predicted cover -------------------------------------------------------------------------------------# 
pred <- readRDS('output/ibm/simulations/ARTR_one_step_ahead_climate_model_cover_per_quad.RDS')

pred$total_size[ is.na(pred$total_size  ) ]  <- 0 

pred <- 
  pred %>% 
  mutate( predicted = rec_area + total_size ) %>%  
  arrange( simulation, Treatment, quad, year)     # prediction for year + 1 

head( obs_cover)
head( pred )

predicted_df  <- merge( obs_cover[, c('year', 'Treatment', 'quad', 'observed', 'Period')], pred[, c('simulation', 'quad', 'Treatment', 'predicted', 'year')] )

head( predicted_df ) 

predicted_pgr <- 
  predicted_df %>% 
  mutate( pgr_predicted = log( predicted)  - log(observed)) %>% 
  group_by( year, Treatment, simulation ) %>% 
  filter( is.finite(pgr_predicted)) %>%
  summarise( mpgr = mean(pgr_predicted)) %>% 
  ungroup(.) %>% 
  group_by(year, Treatment ) %>% 
  summarise( predicted = mean(mpgr) , lcl = quantile ( mpgr, 0.25 ) , ucl = quantile( mpgr , 0.75)) %>% 
  ungroup() %>% 
  mutate( year = year + 1  ) # predictions are for year t + 1 


obs_pgr <- 
  obs_cover %>% 
  ungroup() %>% 
  arrange( Period, Treatment, quad, year ) %>% 
  rename( cover_obs = observed ) %>% 
  mutate( pgr_observed = log(cover_obs) - log(lag(cover_obs,1)), year_diff = year - lag( year , 1 )) %>% 
  filter( year_diff == 1 ) %>% 
  group_by( Period , year , Treatment ) %>% 
  summarise ( observed = mean(pgr_observed))

# merge predicted and observed 

plot_df <- merge(obs_pgr, predicted_pgr)

plot_df_long <- 
  plot_df %>% 
  gather( type, val , observed , predicted ) %>% 
  mutate(lcl = ifelse( type == 'observed' , NA, lcl )) %>% 
  mutate(ucl = ifelse( type == 'observed' , NA, ucl ))
  
library(ggplot2)
ggplot( plot_df_long, aes( x = year, y = val, color = Treatment , shape = type, linetype = type ) ) + geom_point() + geom_line() 

ggplot ( plot_df, aes( x = predicted,  y = observed, color = Period ) ) + 
  geom_point() + 
  geom_smooth(method= 'lm', se = FALSE)

ggplot ( plot_df, aes( x = predicted, y = observed, color = Period ) ) + geom_point() 

ggplot ( plot_df, aes( x = predicted, y = observed, color = Period ) ) + geom_point() 

