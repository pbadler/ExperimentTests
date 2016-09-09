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
library(stringr)

dm <- c(4, 11 ) # Drought months 
im <- c(5, 10 ) # Irrigation months

p.treatments <- c(0.5, 1.5) # Drought and Irrigation adjustments to precip 
t.treatments <- c(1, 1)     # Drought and Irrigation adjustments to temperature 

# make time periods --------------------------------------------------------------------

p1 <- data.frame( Period = 'Modern', year = 2007:2016)
p2 <- data.frame( Period = 'not monitored', year = 1957:2006)
p3 <- data.frame( Period = 'Historical', year = 1926:1956)
periods <- data.frame( rbind( p1, p2, p3 )) 

# ----- read in data --------------------------------------------------------------------#

station_dat <- read.csv('data/USSES_climate_monthly.csv')
season <- read.csv('data/season_table.csv')

# ---process dates----------------------------------------------------------------------#

station_dat$date <-  as.POSIXct( strptime( station_dat$DATE, '%Y%m%d', tz = 'MST')  ) 
station_dat$month <- as.numeric(strftime( station_dat$date, '%m'))
station_dat$year <- as.numeric(strftime( station_dat$date, '%Y'))

# -------- NAs -------------------------------------------------------------------------#

station_dat$TPCP[ station_dat$TPCP == -9999.0 ] <- NA
station_dat$MMNT[ station_dat$MMNT == -9999.0 ] <- NA
station_dat$MMXT[ station_dat$MMXT == -9999.0 ] <- NA
station_dat$MNTM[ station_dat$MNTM == -9999.0 ] <- NA

# set-up aggregate seasonal variables for model ----------------------------------------#

df <- merge( station_dat, season, by = 'month')

df <- 
  df %>% 
  mutate( water_year = year + lag_year ) %>% 
  mutate( quarter = cut(month, 4, labels = paste0('Q', 1:4))) %>%
  select(year, quarter, month, year, season, season_label, precip_seasons, water_year, TPCP, MNTM, MMXT, MMNT)

# ---------- annual average Temperature -------------------------------------------------#
annual_MAT <- 
  df %>% 
  group_by( year ) %>%
  summarise (MAT = mean(MNTM, na.rm = TRUE))

# ---------- annual total precip --------------------------------------------------------#
annual_TPPT <- 
  df %>% 
  group_by( year ) %>% 
  summarise( TPPT = mean(TPCP, na.rm = TRUE), n = n())

# ---------- seasonal average Temperature -----------------------------------------------#

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

# ---------- monthly climate  -----------------------------------------------------------#
# 
# incorporate Drought and Irrigation effects only in specificied months 
# 

monthly_clim <- 
  df %>% 
  mutate(Control = TPCP) %>% 
  select(-TPCP) %>% 
  mutate( Drought    = ifelse( year > 2011 & month %in% c(dm[1]:dm[2]), Control*p.treatments[1], Control) ) %>% 
  mutate( Irrigation = ifelse( year > 2011 & month %in% c(im[1]:im[2]), Control*p.treatments[2], Control) ) %>% 
  gather(Treatment, TPCP, Control, Drought, Irrigation)

# ------------ aggregate monthly climate with Treatment effects by quarter ---------------#

quarterly_clim <-
  monthly_clim %>% 
  gather( var, val ,  MMNT, MMXT, MNTM, TPCP ) %>% 
  group_by(Treatment, var, month) %>% 
  mutate( val = ifelse( year > 1925 & is.na(val), mean(val, na.rm = TRUE), val)) %>% # !!!!!!! fill in missing monthly averages after 1925 with monthly average !!!!!!!! 
  group_by( Treatment, var, year, quarter ) %>%                                      # !!!!!!! note missing values for max, and mean temperature in Feb. 1945.  !!!!!!!!
  summarise( avg = mean(val), ttl = sum(val) ) %>% 
  group_by(var) %>% 
  gather( stat, val, avg, ttl ) %>% 
  filter( (var == 'TPCP' & stat == 'ttl')| (str_detect(pattern = "^M", var) & stat == 'avg')) %>% 
  group_by( Treatment, var, stat) %>% 
  arrange(year, quarter) %>%
  ungroup() %>% 
  unite(var, c(var, stat) , sep = '_')

# -------------------------------------------------------------------------------------------#
# -------------- aggregate monthly climate with Treatment effects by season -----------------#

seasonal_precip <- 
  monthly_clim %>% 
  select ( -MMNT, -MMXT, -MNTM) %>%  
  group_by( Treatment, water_year, precip_seasons ) %>% 
  summarise( TPCP = sum(TPCP) ) %>% 
  group_by( Treatment, precip_seasons ) %>% 
  arrange( water_year, precip_seasons) %>%
  rename( year = water_year ) %>% 
  arrange( year ) %>% 
  rename( l0 = TPCP) %>% 
  mutate( l1 = lag (l0, 1 ) , 
          l2 = lag (l0, 2 ) ) %>% 
  gather( lag, TPCP, l0:l2) %>% 
  ungroup() %>% 
  unite( stat , precip_seasons, lag , sep = '_PRCP_') %>% 
  spread( stat, TPCP )

# -------------- join dfs for variables -------------------------------------------------------#

seasonal_clim <- left_join( seasonal_tmean, seasonal_precip, by = 'year') 

annual_clim <- left_join( annual_TPPT, annual_MAT)

seasonal_clim <- left_join( annual_clim, seasonal_clim)

# -------join periods -------------------------------------------------------------------------#

seasonal_clim <- left_join( seasonal_clim, periods )
monthly_clim <- left_join( df, periods ) 
quarterly_clim <- left_join( quarterly_clim, periods ) 
annual_clim <- left_join( annual_clim, periods ) 

# ------ filter out treatments ----------------------------------------------------------------# 
seasonal_clim <- seasonal_clim %>% 
  filter( !(year < 2012 & Treatment != 'Control')) # remove all Drought and Irrigation treatments prior to 2012 

# quarterly_clim <- quarterly_clim %>% 
#   filter( !(year < 2012 & Treatment != 'Control'))

# -------- output -----------------------------------------------------------------------------#

saveRDS( seasonal_clim, 'data/temp_data/seasonal_climate.RDS')
saveRDS( monthly_clim, 'data/temp_data/monthly_climate.RDS')
saveRDS( quarterly_clim, 'data/temp_data/quarterly_climate.RDS')
saveRDS( annual_clim, 'data/temp_data/annual_climate.RDS')
