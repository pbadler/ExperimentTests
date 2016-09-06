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


dm <- c(4, 11 ) # drought months 
im <- c(5, 10 ) # irrigation months

# make time periods --------------------------------------------------------------------

p1 <- data.frame( period = 'contemporary', year = 2007:2016)
p2 <- data.frame( period = 'not monitored', year = 1956:2006)
p3 <- data.frame( period = 'historical', year = 1926:1955)
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
  df %>% group_by( year ) %>% 
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
# incorporate drought and irrigation effects only in specificied months 
# 

monthly_clim <- 
  df %>% 
  mutate(control = TPCP) %>% 
  select(-TPCP) %>% 
  mutate( drought    = ifelse( month %in% c(dm[1]:dm[2]), control*0.5, control) ) %>% 
  mutate( irrigation = ifelse( month %in% c(im[1]:im[2]), control*1.5, control) ) %>% 
  gather(treatment, TPCP, control, drought, irrigation)

# ------------ aggregate monthly climate with treatment effects by quarter ---------------#

quarterly_climate <-
  monthly_clim %>% 
  gather( var, val ,  MMNT, MMXT, MNTM, TPCP ) %>% 
  group_by( treatment, var, year, quarter ) %>% 
  summarise( avg = mean(val), ttl = sum(val) ) %>% 
  group_by(var) %>% 
  gather( stat, val, avg, ttl ) %>% 
  filter( !(var == 'TPCP' & stat == 'avg' )) %>% 
  group_by( treatment, var, stat) %>% 
  arrange(year, quarter) %>%
  ungroup() %>% 
  unite(var, c(var, stat) , sep = '_')

# -------------------------------------------------------------------------------------------#
# -------------- aggregate monthly climate with treatment effects by season -----------------#

seasonal_precip <- 
  monthly_clim %>% 
  select ( -MMNT, -MMXT, -MNTM) %>%  
  group_by( treatment, water_year, precip_seasons ) %>% 
  summarise( TPCP = sum(TPCP) ) %>% 
  group_by( treatment, precip_seasons ) %>% 
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

# ------ filter out treatments ----------------------------------------------------------------# 
seasonal_clim <- seasonal_clim %>% 
  filter( !(year < 2012 & treatment != 'control')) # remove all drought and irrigation treatments prior to 2012 

# -------- output -----------------------------------------------------------------------------#

saveRDS( seasonal_clim, 'data/temp_data/seasonal_climate.RDS')
saveRDS( df, 'data/temp_data/monthly_climate.RDS')
saveRDS( quarterly_climate, 'data/temp_data/quarterly_climate.RDS')
