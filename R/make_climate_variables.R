rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(lme4)
library(zoo)

station_dat <- read.csv('data/USSES_climate.csv')

season <- readRDS('data/temp_data/season.RDS')

# ---------------------------------------------------------------------------------------
station_dat$date <-  as.POSIXct( strptime( station_dat$DATE, '%Y%m%d', tz = 'MST')  ) 
station_dat$doy <- as.numeric(strftime( station_dat$date, '%j'))
station_dat$month <- as.numeric(strftime( station_dat$date, '%m'))
station_dat$year <- as.numeric(strftime( station_dat$date, '%Y'))

station_dat$PRCP[ station_dat$PRCP == -9999.0 ] <- NA
station_dat$TMAX[ station_dat$TMAX == -9999.0 ] <- NA
station_dat$TMIN[ station_dat$TMIN == -9999.0 ] <- NA

# set-up aggregate seasonal variables for model ----------------------------------------- 

df <- merge( station_dat, season, by = 'month')

df <- 
  df %>% 
  mutate( TMEAN = ( TMAX + TMIN ) / 2 ) %>% 
  mutate( water_year = year + lag_year ) %>% 
  select(date, doy, month, year, season, season_label, precip_seasons, water_year, PRCP, TMEAN)

seasonal_tmean <- 
  df %>% 
  group_by(year, season_label) %>% 
  summarise( l0 = mean(TMEAN, na.rm = TRUE) )

seasonal_tmean <- 
  seasonal_tmean %>% 
  ungroup() %>% 
  group_by(season_label) %>% 
  arrange( season_label, year ) %>% 
  mutate( l1 = lag ( l0, 1 ), 
          l2 = lag ( l0, 2 ) ) %>% 
  gather( lag, TMEAN, l0:l2 ) %>% 
  ungroup( ) %>% 
  unite( stat,  season_label, lag , sep = '_TMEAN_') %>% 
  spread( stat, TMEAN )

monthly_precip <- 
  df %>% 
  group_by(year, water_year, precip_seasons, month) %>% 
  summarise(control = sum(PRCP, na.rm = TRUE)) %>% 
  mutate( drought    = ifelse( month %in% c(3:11), control*0.5, control) ) %>% 
  mutate( irrigation = ifelse( month %in% c(4:10), control*1.5, control) ) %>%
  gather( treatment, )
  
  


seasonal_precip <- 
  df %>% 
  group_by(water_year, precip_seasons ) %>% 
  summarise(l0 = sum(PRCP)) %>% 
  rename( year = water_year ) %>% 
  group_by( precip_seasons ) %>% 
  arrange( precip_seasons, year ) %>% 
  mutate( l1 = lag (l0, 1 ) , 
          l2 = lag (l0, 2 ) ) %>% 
  gather( lag, PRCP, l0:l2) %>% 
  ungroup() %>% 
  unite( stat , precip_seasons, lag , sep = '_PRCP_') %>% 
  spread( stat, PRCP )


seasonal_clim <- left_join( seasonal_tmean, seasonal_precip, by = 'year') 

saveRDS( seasonal_clim, 'data/temp_data/seasonal_climate.RDS')

