#######################################################################################
#
# Setup seasonal climate variables for demographic rate models 
#
#######################################################################################

rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(lme4)
library(zoo)

station_dat <- read.csv('data/USSES_climate.csv')
station_dat <- read.csv('data/USSES_climate_monthly.csv')

season <- read.csv('data/season_table.csv')

# ---------------------------------------------------------------------------------------
station_dat$date <-  as.POSIXct( strptime( station_dat$DATE, '%Y%m%d', tz = 'MST')  ) 
#station_dat$doy <- as.numeric(strftime( station_dat$date, '%j'))
station_dat$month <- as.numeric(strftime( station_dat$date, '%m'))
station_dat$year <- as.numeric(strftime( station_dat$date, '%Y'))

station_dat$TPCP[ station_dat$TPCP == -9999.0 ] <- NA
station_dat$MMNT[ station_dat$MMNT == -9999.0 ] <- NA
station_dat$MMXT[ station_dat$MMXT == -9999.0 ] <- NA
station_dat$MNTM[ station_dat$MNTM == -9999.0 ] <- NA

# set-up aggregate seasonal variables for model ----------------------------------------- 

df <- merge( station_dat, season, by = 'month')

df <- 
  df %>% 
  mutate( water_year = year + lag_year ) %>% 
  select(year, month, year, season, season_label, precip_seasons, water_year, TPCP, MNTM, MMXT, MMNT)

annual_MAT <- 
  df %>% 
  group_by( year ) %>%
  summarise (MAT = mean(MNTM, na.rm = TRUE))

annual_TPPT <- 
  df %>% group_by( year ) %>% 
  summarise( TPPT = mean(TPCP, na.rm = TRUE), n = n())

seasonal_tmean <- 
  df %>% 
  mutate(year = ifelse(month == 12 , year + 1, year  )) %>% # account for December 
  group_by(year, season_label) %>% 
  summarise( l0 = mean(MNTM, na.rm = TRUE) )

seasonal_tmean <- 
  seasonal_tmean %>% 
  ungroup() %>% 
  group_by(season_label)%>% 
  arrange( season_label, year ) %>% 
  mutate( l1 = lag ( l0, 1 ), 
          l2 = lag ( l0, 2 ) ) %>% 
  gather( lag, MNTM, l0:l2 ) %>% 
  ungroup( ) %>% 
  unite( stat,  season_label, lag , sep = '_TMEAN_') %>% 
  spread( stat, MNTM )

monthly_precip <- 
  df %>% 
  group_by(year, water_year, precip_seasons, month) %>% 
  summarise(control = sum(TPCP, na.rm = TRUE)) %>% 
  mutate( drought    = ifelse( month %in% c(3:11), control*0.5, control) ) %>% 
  mutate( irrigation = ifelse( month %in% c(4:10), control*1.5, control) ) %>% 
  gather(treatment, TPCP, control, drought, irrigation)
  
seasonal_precip <- 
  df %>% 
  group_by(year, water_year, precip_seasons, month) %>% 
  summarise(Control = sum(TPCP, na.rm = TRUE)) %>% 
  mutate( Drought    = ifelse( month %in% c(3:11), Control*0.5, Control) ) %>% 
  mutate( Irrigation = ifelse( month %in% c(4:10), Control*1.5, Control) ) %>% 
  gather(Treatment, TPCP, Control, Drought, Irrigation) %>%
  group_by(water_year, precip_seasons , Treatment ) %>% 
  summarise(l0  = sum(TPCP, na.rm = TRUE)) %>% 
  rename( year = water_year ) %>% 
  group_by( Treatment, precip_seasons ) %>% 
  arrange( year ) %>% 
  mutate( l1 = lag (l0, 1 ) , 
          l2 = lag (l0, 2 ) ) %>% 
  gather( lag, TPCP, l0:l2) %>% 
  ungroup() %>% 
  unite( stat , precip_seasons, lag , sep = '_PRCP_') %>% 
  spread( stat, TPCP )

seasonal_clim <- left_join( seasonal_tmean, seasonal_precip, by = 'year') 

annual_clim <- left_join( annual_TPPT, annual_MAT)

seasonal_clim <- left_join( annual_clim, seasonal_clim)

# join periods

p1 <- data.frame( period = 'contemporary', year = 2007:2016)
p2 <- data.frame( period = 'not monitored', year = 1956:2006)
p3 <- data.frame( period = 'historical', year = 1926:1955)

seasonal_clim <- left_join( seasonal_clim, data.frame( rbind(p1, p2, p3)) )

seasonal_clim <- seasonal_clim %>% filter( !(year < 2012 & Treatment != 'Control'))

saveRDS( seasonal_clim, 'data/temp_data/seasonal_climate.RDS')
saveRDS( df, 'data/temp_data/monthly_climate.RDS')
