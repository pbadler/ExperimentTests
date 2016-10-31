# calculate treatment effects compared to control soil moisture 

rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(lme4)
library(zoo)

myVWC <- readRDS('~/driversdata/data/idaho_modern/soil_moisture_data/data/processed_data/decagon_data_with_station_data.RDS')
#spotVWC <- readRDS('~/driversdata/data/idaho_modern/soil_moisture_data/data/processed_data/spring_spot_measurements.RDS')
swVWC <- read.csv('data/temp_data/daily_VWC.csv')

# --------------------------------------------------------------------------------------#

myVWC <- myVWC %>% 
  mutate( v = ifelse(measure == 'VWC', v*100, v)) # convert to percent 

myVWC$year <- as.numeric( strftime( as.Date( myVWC$simple_date), '%Y' ) )

df_soil_moist <- 
  myVWC %>% 
  filter( bad_values == 0 , stat == 'raw', measure == 'VWC') %>% 
  filter( !is.na(rainfall)) %>% 
  filter( !(plot == 16) ) %>%           ##### Drop plot 16  
  group_by( year, PrecipGroup, Treatment, season, simple_date, rainfall, depth ) %>% 
  summarise( avg_VWC = mean(v, na.rm = TRUE))  %>% 
  group_by(PrecipGroup, depth) %>% 
  mutate( avg_VWC = scale(avg_VWC)) %>%  # scale within Precip Group and Depth 
  spread( Treatment, avg_VWC) %>%
  mutate( Drought = Drought - Control, Irrigation = Irrigation - Control )

# ---process dates----------------------------------------------------------------------#

swVWC$date <- as.POSIXct(strptime( paste( swVWC$Year, swVWC$DOY, sep = '-') , '%Y-%j'))
swVWC$year <- swVWC$Year
swVWC$month <- as.numeric(strftime( swVWC$date, '%m'))

# set-up aggregate seasonal variables for model ----------------------------------------#

swVWC <- 
  swVWC %>% 
  gather(layer, VWC, Lyr_1:Lyr_6) 

# ----------- aggregate layers to two layers, shallow and deep  --------------------# 

swVWC$layer [ grep( '[1-3]+' , swVWC$layer)  ] <- 1 
swVWC$layer [ grep( '[4-6]+' , swVWC$layer ) ] <- 2

swVWC <- swVWC %>% 
  group_by( date, layer ) %>% 
  summarise( modelVWC = mean(VWC)*100 )

# ----------------- split into depth layers -----------------------------------------# 

VWC_df <- split(df_soil_moist, df_soil_moist$depth)

# model the treatment effect within each group by the control VWC -------------------# 

# fit models 
m5cmDrought <- lm( data = VWC_df$`5 cm deep`, Drought ~ rainfall*season)
m5cmIrrigation <- lm(data = VWC_df$`5 cm deep`, Irrigation ~ rainfall*season)

m25cmDrought <- lm( data = VWC_df$`25 cm deep`, Drought ~ rainfall*season)
m25cmIrrigation <- lm(data = VWC_df$`25 cm deep`, Irrigation ~  rainfall*season)

# select models 

library(MASS)
m5cmDrought <- stepAIC(m5cmDrought,scope=list(upper=~.,lower=~1),trace=T)
m5cmIrrigation <- stepAIC(m5cmIrrigation,scope=list(upper=~.,lower=~1),trace=T)

m25cmDrought <- stepAIC(m25cmDrought,scope=list(upper=~.,lower=~1),trace=T)
m25cmIrrigation <- stepAIC(m25cmIrrigation,scope=list(upper=~.,lower=~1),trace=T)

#

daily_avg_VWC <- 
  do.call(rbind, VWC_df) %>% 
  gather( Treatment, VWC, Control:Irrigation) %>% 
  group_by(rainfall, year, season, simple_date, depth,Treatment) %>% 
  summarise( avg_VWC = mean(VWC, na.rm = TRUE) ) %>% 
  spread( Treatment, avg_VWC ) 

daily_avg_VWC <- split(daily_avg_VWC, factor( daily_avg_VWC$depth))

daily_avg_VWC$`5 cm deep`$Dpred <- predict( m5cmDrought, daily_avg_VWC$`5 cm deep`)
daily_avg_VWC$`5 cm deep`$Ipred <- predict( m5cmIrrigation, daily_avg_VWC$`5 cm deep`)

daily_avg_VWC$`25 cm deep`$Dpred <- predict( m25cmDrought, daily_avg_VWC$`25 cm deep`)
daily_avg_VWC$`25 cm deep`$Ipred <- predict( m25cmIrrigation, daily_avg_VWC$`25 cm deep`)

ggplot( daily_avg_VWC$`5 cm deep`, aes( x = simple_date, y = Control)) +
  geom_line() +
  geom_line(aes( y = Irrigation + Control) , color = 'blue', alpha = 0.2) + 
  geom_line(aes( y = Ipred + Control), color = 'blue', linetype = 2) + 
  geom_line(aes( y = Drought + Control) , color = 'red', alpha = 0.2) + 
  geom_line(aes( y = Dpred + Control), color = 'red', linetype = 2) 


# predict the treatment effects from the scaled model VWC ----------------------------# 
daily_clim <- readRDS('~/driversdata/data/idaho_modern/soil_moisture_data/data/processed_data/daily_station_dat_rainfall.RDS')
seasons <- read.csv('~/driversdata/data/idaho_modern/soil_moisture_data/data/season_table.csv')

daily_clim$date <- as.Date(daily_clim$date ) 
daily_clim$month <- as.numeric( strftime(daily_clim$date, '%m') ) 
swVWC$date <- as.Date( swVWC$date)

daily_clim <- left_join(daily_clim, seasons, by = 'month')
swVWC <- left_join( swVWC, daily_clim , by = 'date') 

swVWC <- split( swVWC , swVWC$layer ) 

swVWC$`1`$Control <- scale( swVWC$`1`$modelVWC ) # standardize control SWC 
swVWC$`2`$Control <- scale( swVWC$`2`$modelVWC ) 

swVWC$`1`$Drought <- predict(m5cmDrought,  swVWC$`1`)
swVWC$`1`$Irrigation <- predict(m5cmIrrigation, swVWC$`1`)

swVWC$`2`$Drought <- predict(m25cmDrought,  swVWC$`2`)
swVWC$`2`$Irrigation <- predict(m25cmIrrigation, swVWC$`2`)

# unscale the Control Drought and Irrigation VWC --------------------------------------------------------  #

swVWC$`1`$center <- attr(swVWC$`1`$Control, 'scaled:center')
swVWC$`1`$scale  <- attr(swVWC$`1`$Control, 'scaled:scale')

swVWC$`2`$center <- attr(swVWC$`2`$Control, 'scaled:center')
swVWC$`2`$scale  <- attr(swVWC$`2`$Control, 'scaled:scale')

swVWC <- do.call(rbind, swVWC)

swVWC <- 
  swVWC %>% 
  mutate( Drought = Control + Drought, Irrigation = Control + Irrigation) %>% 
  gather( Treatment, VWC, Control:Irrigation) 

swVWC$VWC_raw <- swVWC$VWC*swVWC$scale + swVWC$center

my_colors <- c('#66c2a5','#fc8d62','#8da0cb')

ggplot( swVWC, aes( x = date, y = VWC_raw, color = Treatment)) + 
  geom_line() + 
  scale_color_manual(values = my_colors) + 
  xlim( as.Date( c('2014-01-01', '2015-01-01'))) + 
  facet_wrap( ~ layer )

ggplot( swVWC, aes( x = date, y = VWC, color = Treatment)) + 
  geom_line() + 
  scale_color_manual(values = my_colors) + 
  xlim( as.Date( c('2014-01-01', '2015-01-01'))) + 
  facet_wrap( ~ layer )

saveRDS(swVWC, 'data/temp_data/daily_swVWC_treatments.RDS')

