rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(lme4)
library(zoo)

df <- readRDS(file = 'data/temp_data/decagon_data_corrected_values.RDS')
station_dat <- read.csv('data/USSES_climate.csv')

# ---------------------------------------------------------------------------------------

df <- df %>% filter( stat == 'raw', bad_values == 0 ) %>% mutate(v = ifelse(measure == 'VWC' , v*100, v))

df$depth <- factor( df$depth , levels = c('air temperature', '5 cm deep', '25 cm deep') , order = TRUE ) 
df$Treatment <- factor(df$Treatment, levels = c('Drought', 'Control', 'Irrigation'), order = TRUE)

df$datetime <- df$date

station_dat$date <-  as.POSIXct( strptime( station_dat$DATE, '%Y%m%d')  ) 
station_dat$date <- as.Date(station_dat$date)

station_dat <- station_dat %>% 
  mutate( TMEAN = ( TMAX + TMIN ) / 2 ) %>% 
  filter(date > '2011-01-01') %>% 
  select(date, PRCP, TMEAN)

station_dat$PRCP[ station_dat$PRCP == -9999.0 ] <- NA
station_dat$TMEAN[ station_dat$TMEAN == -9999.0 ] <- NA

station_dat <- station_dat %>% 
  mutate( rainfall = rollapply(PRCP, 2, sum, fill = 0, na.rm = TRUE, align = 'right') ) %>%
  mutate( rainfall = ifelse( rainfall > 1.0 & TMEAN > 3 & !is.na(rainfall), 'rainy', 'not rainy'))

df$month <- strftime(df$datetime, format = '%m' )
df$month <- as.numeric( df$month)
df$hour <- strftime( df$datetime, format = '%H')
df$hour <- as.numeric( df$hour)
df$date <- as.Date(df$datetime)

season <- data.frame ( month = 1:12, season = c('winter', 'winter', 'spring', 'spring', 'spring', 'summer', 'summer', 'summer', 'fall', 'fall', 'fall', 'winter'))

season$season <- factor( season$season, levels = c('spring', 'summer', 'fall', 'winter'), order = TRUE)

tod <- data.frame( hour = 1:24, tod = cut(1:24, breaks = c(0, 6, 19, 24)) )
levels(tod$tod ) <- c('night', 'day', 'night')

df <- merge( df, season, by = 'month')
df <- merge( df, tod, by = 'hour')

df <- df %>% left_join( station_dat, by = 'date') 

# summarize treatment differences:  -----------------------------------------------------------------------------------

plot_vals <- df %>% 
  filter( !is.na(v), bad_values == 0 ) %>% 
  group_by(Treatment, season, depth, measure ) %>% 
  summarise( avg = mean(v), stddev = sd(v), max = max(v), min = min(v), n = n(), ci = 1.96*(stddev/sqrt(n)), uci = avg + ci, lci = avg - ci ) 

plot_vals_rainy <- df %>% 
  filter( !is.na(v), bad_values == 0, month > 3, month < 12, measure == 'VWC' , !is.na(rainfall)) %>% 
  group_by(Treatment, season, depth, measure, rainfall ) %>% 
  summarise( avg = mean(v), stddev = sd(v), max = max(v), min = min(v), n = n(), ci = 1.96*(stddev/sqrt(n)), uci = avg + ci, lci = avg - ci )

group_vals <- df %>% 
  filter( !is.na(v), bad_values == 0 ) %>% 
  group_by(PrecipGroup, Treatment, season, depth, measure) %>% 
  summarise( avg = mean(v), stddev = sd(v), max = max(v), min = min(v), n = n(), ci = 1.96*(stddev/sqrt(n)), uci = avg + ci, lci = avg - ci ) 

brplt <- ggplot(plot_vals, aes( x = season, y = avg, fill = Treatment)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar( aes(ymax = uci, ymin = lci), position = 'dodge') + 
  facet_wrap( ~ depth ) 

brplt_rainy <- brplt + facet_wrap(  depth ~ rainfall) 

plot_T_vals <- df %>% 
  filter( !is.na(v), bad_values == 0, measure == 'C', depth != '25 cm deep') %>%
  group_by(Treatment, season, tod, depth) %>% 
  summarise( avg = mean(v), stddev = sd(v), max = max(v), min = min(v), n = n(), ci = 1.96*(stddev/sqrt(n)), uci = avg + ci, lci = avg - ci ) 

brplt_tod <- ggplot(plot_vals, aes( x = season, y = avg, fill = Treatment)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar( aes(ymax = uci, ymin = lci), position = 'dodge') + 
  facet_wrap( depth ~  tod ) 

plot_groups <- ggplot( group_vals, aes ( x = factor( PrecipGroup ) , y = avg, fill = Treatment )) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar( aes(ymax = uci, ymin = lci), position = 'dodge') + 
  facet_wrap( ~ depth ) 

brplt %+% subset(plot_vals, measure == 'VWC' ) + ylab( 'Soil Moisture')

brplt_rainy %+% subset(plot_vals_rainy , measure == 'VWC') + ylab( 'Volumetric Soil Moisture (%)')

brplt %+% subset(plot_vals, measure == 'C' & depth != '25 cm deep') + ylab( 'Temperature (C)')

brplt_tod %+% plot_T_vals + ylab( 'Temperature (C)')

plot_groups %+% subset( group_vals, measure == 'VWC' & season == 'spring') + ylab ( "average spring soil moisture") + xlab ( "plot group")
  
# model formula and models ----------------------------------------------------------------------- 

df <- df %>% 
  mutate ( unique_port = paste0( plot, '.', port))

df_5cm_soil <- df %>% 
  filter( bad_values == 0 , stat == 'raw', measure == 'VWC', depth == '5 cm deep')

df_25cm_soil <- df %>% 
  filter( bad_values == 0, stat == 'raw', measure == 'VWC', depth == '25 cm deep')

df_air_temp <- df %>% 
  filter( bad_values == 0 , stat == 'raw', measure == 'C', depth == 'air temperature') 

df_soil_temp <- df %>% 
  filter( bad_values == 0 , stat == 'raw', measure == 'C', depth == '5 cm deep')

basic_form <- formula(v ~ (1|unique_port) + (1|datetime) + factor( PrecipGroup )  + Treatment*season )
VWC_form <- formula( v ~ (1|unique_port) + (1|datetime) + factor( PrecipGroup) + Treatment*rainfall ) 

m_air <- lmer(data =  df_air_temp, formula = basic_form)
m_soil <- lmer(data = df_soil_temp, formula = basic_form)

m_5cm <- lmer( data = df_5cm_soil, formula = update( basic_form, . ~ . + (1|plot) + Treatment*season*rainfall))
m_25cm <- lmer( data = df_25cm_soil, formula = update( basic_form, . ~ . + (1|plot) + Treatment*season*rainfall)) 

m_5cm_spring <- lmer( data = subset(df_5cm_soil, season == 'spring'), formula = VWC_form)
m_5cm_summer <- lmer( data = subset(df_5cm_soil, season == 'summer'), formula = VWC_form ) 

summary(m_5cm_spring)
summary(m_5cm_summer)

summary(m_air)
summary(m_soil)
summary(m_5cm) 
summary(m_25cm)


# ------------------------------------------------------------------------------------------------- 

# make prediction df to view lmer effects 

pred_df <- expand.grid( PrecipGroup = c(1,3,4,6), Treatment = unique( df$Treatment), season = levels(df$season), rainfall = levels( factor(df$rainfall) ) )  

pred_df$air_temp_pred <- predict(m_air, pred_df, re.form = NA)
pred_df$VWC_5cm_pred <- predict( m_5cm, pred_df, re.form = NA )
pred_df$VWC_25cm_pred <- predict( m_25cm, pred_df, re.form = NA)
pred_df$VWC_5_cm_spring <- predict( m_5cm_spring, pred_df, re.form = NA)

pred_stats_T <- pred_df %>% 
  group_by( Treatment, season) %>% 
  summarise( air_pred = mean(air_temp_pred)) 

pred_stats_VWC <- pred_df %>% 
  group_by(season, rainfall, Treatment ) %>% 
  summarise( VWC_5cm_pred = mean(VWC_5cm_pred), VWC_25cm_pred = mean(VWC_25cm_pred)) 

pred_plot <- 
  ggplot( pred_stats_T, aes( x = season, y = air_pred, fill = Treatment )) + 
  geom_bar( stat = 'identity', position = 'dodge') + 
  ggtitle( 'Treatment effects predicted by lmer model') + 
  ylab( 'Air temperature (C)')

pred_plot_VWC_5 <-
  ggplot( pred_stats_VWC, aes( x = season, y = VWC_5cm_pred, fill = Treatment )) + 
  geom_bar( stat = 'identity', position = 'dodge') + 
  facet_wrap( ~ rainfall ) + 
  ggtitle( 'Treatment effects predicted by lmer model') + 
  ylab( '5 cm soil moisture (%)')
  
pred_plot_VWC_25 <- 
  ggplot( pred_stats_VWC, aes ( x = season, y = VWC_25cm_pred, fill = Treatment)) + 
  geom_bar( stat = 'identity', position = 'dodge') + 
  facet_wrap( ~ rainfall ) + 
  ggtitle( 'Treatment effects predicted by lmer model') + 
  ylab( '25 cm soil moisture (%)')

# print plots in one pdf ---------------------------------------------------------------------- 

pdf( "figures/summary_of_treatment_effects.pdf" , height = 8, width = 10 ) 

print( brplt %+% subset(plot_vals, measure == 'VWC' ) + ylab( 'Volumetric Soil Moisture (%)') ) 

print( brplt %+% subset(plot_vals, measure == 'C' & depth != '25 cm deep') + ylab( 'Temperature (C)') ) 

print( brplt_rainy %+% subset(plot_vals_rainy , measure == 'VWC') + ylab( 'Volumetric Soil Moisture (%)') )

print( brplt_tod %+% plot_T_vals + ylab( 'Temperature (C)') ) 

print ( plot_groups %+% subset( group_vals, measure == 'VWC' & season == 'spring') + ylab ( "Spring soil moisture (%)") + xlab ( "plot group") ) 

print( pred_plot ) 

print( pred_plot_VWC_5)

print( pred_plot_VWC_25)

dev.off() 


