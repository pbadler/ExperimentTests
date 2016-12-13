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

df$true_cov <- exp(df$logarea.t0)

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
  arrange( simulation, Treatment, quad, year) %>% 
  mutate( year = year + 1 )                        #### plus one to get line up prediction with observed year 

# %>% 
#   mutate( change = log( cover ) - log(lag(cover, 1) ), year_diff = year - lag(year, 1)) %>% 
#   group_by(quad, year, Treatment) %>%
#   filter( year_diff == 1) %>% 
#   filter( is.finite( change )) %>%
#   summarise( cover_pred = mean(cover),  predicted = mean(change), ucl = quantile( change, 0.75), lcl = quantile( change, 0.25) )

head( obs_cover)
head( pred )

predicted_df  <- merge( obs_cover[, c('year', 'Treatment', 'quad', 'observed', 'Period')], pred[, c('simulation', 'quad', 'Treatment', 'predicted', 'year')] )

predicted_df %>% 
  mutate( pgr_predicted =  predicted - lag (observed, 1)) %>% 
  group_by( year, Treatment , quad ) %>% 
  filter( )
  summarise( mean(pgr_predicted), ucl = quantile( pgr_predicted, 0.75), lcl = quantile( pgr_predicted, 0.25 )) 

obs_change <- 
  obs_cover %>% 
  ungroup() %>% 
  arrange( Treatment, quad, year ) %>% 
  rename( cover_obs = observed ) %>% 
  mutate( observed = log(cover_obs) - log(lag(cover_obs,1)), year_diff = year - lag( year , 1 )) %>% 
  filter( year_diff == 1 ) 


# merge predicted and observed 

plot_df <- merge( obs_change[ , c('Period', 'Group', 'Treatment', 'quad', 'cover_obs', 'observed', 'year')], pred, all.x = T)

ggplot ( plot_df, aes( x = cover_pred, y = cover_obs, color = Period ) ) + geom_point() 

ggplot ( plot_df, aes( x = predicted, y = observed, color = Period ) ) + geom_point() 

plot_df 
