# calculate treatment effects compared to control soil moisture 

rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(lme4)
library(zoo)

myVWC <- readRDS('~/driversdata/data/idaho_modern/soil_moisture_data/data/processed_data/decagon_data_with_station_data.RDS')
swVWC <- read.csv('data/temp_data/daily_VWC.csv')
daily_clim <- readRDS('~/driversdata/data/idaho_modern/soil_moisture_data/data/processed_data/daily_station_dat_rainfall.RDS')
seasons <- read.csv('~/driversdata/data/idaho_modern/soil_moisture_data/data/season_table.csv')

# --------------------------------------------------------------------------------------#
myVWC <- myVWC %>% 
  mutate( v = ifelse(measure == 'VWC', v*100, v)) # convert to percent 

myVWC$year <- as.numeric( strftime( as.Date( myVWC$simple_date), '%Y' ) )

# Steps for treatment effects standardization: 
  # 1. aggregate soil moisture by PrecipGroup, Treatment, depth and date 
  # 2. then standardize soil moisture measurements within PrecipGroup and depth
  # 3. then subtract standardized Drought and Irrigation from Control to get 
  #       standardized treatment effects. 
  # 4. output 

VWC_weights <- 
  myVWC %>% 
  filter( bad_values == 0 , stat == 'raw', measure == 'VWC') %>% 
  filter( !is.na(rainfall)) %>% 
  filter( !(plot == 16 ) ) %>%           ##### Drop plot 16  
  group_by( year, PrecipGroup, simple_date) %>%  
  summarise( weight = n_distinct(unique_position))

df_soil_moist <- 
  myVWC %>% 
  filter( bad_values == 0 , stat == 'raw', measure == 'VWC') %>% 
  filter( !is.na(rainfall)) %>% 
  filter( !(plot == 16 ) ) %>%           ##### Drop plot 16  
  group_by( year, PrecipGroup, Treatment, season, simple_date, rainfall) %>%  
  summarise( avg_VWC = mean(v, na.rm = TRUE)) %>% 
  group_by(PrecipGroup) %>% 
  mutate( avg_VWC = scale(avg_VWC)) %>%  # scale within Precip Group and Depth 
  spread( Treatment, avg_VWC) %>%
  mutate( Drought = Drought - Control, Irrigation = Irrigation - Control ) %>% 
  arrange( PrecipGroup, simple_date)

df_soil_moist <- merge( df_soil_moist, VWC_weights)

# ---process dates----------------------------------------------------------------------#

swVWC$date <- as.POSIXct(strptime( paste( swVWC$Year, swVWC$DOY, sep = '-') , '%Y-%j'))
swVWC$year <- swVWC$Year
swVWC$month <- as.numeric(strftime( swVWC$date, '%m'))

# set-up aggregate seasonal variables for model ----------------------------------------#

swVWC <- 
  swVWC %>% 
  gather(layer, VWC, Lyr_1:Lyr_6) %>% 
  filter( layer %in% c('Lyr_1', 'Lyr_2', 'Lyr_3'))

# ----------- aggregate layers to two layers, shallow and deep  --------------------# 
swVWC <- swVWC %>% 
  group_by( date) %>% 
  summarise( modelVWC = mean(VWC)*100 )

# ----------------- split into depth layers -----------------------------------------# 
df_soil_moist$type <- 'logger'
VWC_df <- df_soil_moist

spotVWC <- readRDS('data/temp_data/spotVWC.RDS')
spotVWC$type <- 'spot'
spotVWC$year <- strftime( spotVWC$date, '%Y')

VWC_df <- rbind( VWC_df, spotVWC[, -1])

# 
#  Model the standardized treatment differences. 
#  Treatment effects vary with season and whether it's a rainy period. 
#  See definition of rainy periods in the datadrivers data preparation scripts. 
# 

# fit models 
mDrought <- lm( data = VWC_df, Drought ~ rainfall*season , weights = VWC_df$weight)
mIrrigation <- lm(data = VWC_df, Irrigation ~ rainfall*season , weights = VWC_df$weight)

# select models 

library(MASS)
mDrought <- stepAIC(mDrought,scope=list(upper=~.,lower=~1),trace=T)
mIrrigation <- stepAIC(mIrrigation,scope=list(upper=~.,lower=~1),trace=T)

mDrought <- lmer( update( formula( mDrought), . ~ . + (1|simple_date) + (1|PrecipGroup)), data = VWC_df, weights = VWC_df$weight)
mIrrigation <- lmer( update( formula( mIrrigation), . ~ . + (1|simple_date) + (1|PrecipGroup)), data = VWC_df, weights = VWC_df$weight)

#
daily_avg_VWC <- 
  VWC_df %>% 
  gather( Treatment, VWC, Control:Irrigation) %>% 
  group_by(rainfall, year, season, simple_date,Treatment) %>% 
  summarise( avg_VWC = mean(VWC, na.rm = TRUE) ) %>% 
  spread( Treatment, avg_VWC ) 

daily_avg_VWC$Dpred <- predict( mDrought, daily_avg_VWC, re.form = NA)
daily_avg_VWC$Ipred <- predict( mIrrigation, daily_avg_VWC, re.form = NA)

# daily_avg_VWC$`25 cm deep`$Dpred <- predict( m25cmDrought, daily_avg_VWC$`25 cm deep`)
# daily_avg_VWC$`25 cm deep`$Ipred <- predict( m25cmIrrigation, daily_avg_VWC$`25 cm deep`)

pdf( 'figures/treatment_effects_modeled_vs_observed_soil_moisture.pdf', width = 11, height = 5)
print( 
    ggplot( daily_avg_VWC, aes( x = simple_date, y = Control)) +
      geom_line() +
      geom_line(aes( y = Irrigation + Control) , color = 'blue', alpha = 0.2) + 
      geom_line(aes( y = Ipred + Control), color = 'blue', linetype = 2) + 
      geom_line(aes( y = Drought + Control) , color = 'red', alpha = 0.2) + 
      geom_line(aes( y = Dpred + Control), color = 'red', linetype = 2) 
)
dev.off()

# predict the treatment effects from the scaled model VWC ----------------------------# 
daily_clim$date <- as.Date(daily_clim$date ) 
swVWC$date <- as.Date( swVWC$date)
swVWC$month <- as.numeric( strftime(swVWC$date, '%m') )
swVWC$year <- as.numeric( strftime(swVWC$date, '%Y'))

daily_clim <- 
  daily_clim %>% 
  ungroup() %>%
  dplyr::select( -year)

swVWC <- left_join(swVWC, seasons, by = 'month')
swVWC <- left_join( swVWC, daily_clim, by = c('date')) 

swVWC$Control <- scale( swVWC$modelVWC ) # standardize control SWC 
swVWC$Drought <- predict(mDrought,  swVWC, re.form = NA)
swVWC$Irrigation <- predict(mIrrigation, swVWC, re.form = NA)

# unscale the Control Drought and Irrigation VWC --------------------------------------------------------  #

swVWC$center <- attr(swVWC$Control, 'scaled:center')
swVWC$scale  <- attr(swVWC$Control, 'scaled:scale')

swVWC <- 
  swVWC %>% 
  mutate( Drought = Control + Drought, Irrigation = Control + Irrigation) %>% 
  gather( Treatment, VWC, Control:Irrigation) 

swVWC$VWC_raw <- swVWC$VWC*swVWC$scale + swVWC$center

my_colors <- c('#66c2a5','#fc8d62','#8da0cb')

pdf( 'figures/modeled_soilwat_soil_moisture_example.pdf', width = 8, height = 6)
print( 
  ggplot( swVWC, aes( x = date, y = VWC_raw, color = Treatment)) + 
    geom_line() + 
    scale_color_manual(values = my_colors) + 
    xlim( as.Date( c('2016-01-01', '2016-10-01')))
)
dev.off()
print( 
  ggplot( swVWC, aes( x = date, y = VWC, color = Treatment)) + 
    geom_line() + 
    scale_color_manual(values = my_colors) + 
    xlim( as.Date( c('2014-01-01', '2015-01-01'))) 
)
saveRDS(swVWC, 'data/temp_data/daily_swVWC_treatments.RDS')

