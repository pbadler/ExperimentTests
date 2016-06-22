rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(lme4)
library(zoo)

df <- readRDS(file = 'data/temp_data/decagon_data_corrected_values.RDS')

# -------------------------------------------------------------------------------------------------------------- 

df$month <- strftime(df$date, format = '%m' )
df$month <- as.numeric( df$month)
df$hour <- strftime( df$date, format = '%H')
df$hour <- as.numeric( df$hour)
df$day <- as.Date(df$date)

season <- data.frame ( month = 1:12, season = c('winter', 'winter', 'spring', 'spring', 'spring', 'summer', 'summer', 'summer', 'fall', 'fall', 'fall', 'winter'))
tod <- data.frame( hour = 1:24, tod = cut(1:24, breaks = c(0, 6, 19, 24)) )
levels(tod$tod ) <- c('night', 'day', 'night')

df <- merge( df, season, by = 'month')
df <- merge( df, tod, by = 'hour')

# summarize treatment differences: 

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

plot_vals <- df %>% 
  filter( !is.na(v), bad_values == 0 ) %>% 
  group_by(Treatment, season, depth, measure ) %>% 
  summarise( avg = mean(v), stddev = sd(v), max = max(v), min = min(v), n = n(), ci = 1.96*(stddev/sqrt(n)), uci = avg + ci, lci = avg - ci ) 

group_vals <- df %>% 
  filter( !is.na(v), bad_values == 0 ) %>% 
  group_by(PrecipGroup, Treatment, season, depth, measure) %>% 
  summarise( avg = mean(v), stddev = sd(v), max = max(v), min = min(v), n = n(), ci = 1.96*(stddev/sqrt(n)), uci = avg + ci, lci = avg - ci ) 

brplt <- ggplot(plot_vals, aes( x = season, y = avg, fill = Treatment)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_errorbar( aes(ymax = uci, ymin = lci), position = 'dodge') + 
  facet_wrap( ~ depth ) 

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

brplt %+% subset(plot_vals, measure == 'C' & depth != '25 cm deep') + ylab( 'Temperature (C)')

brplt_tod %+% plot_T_vals + ylab( 'Temperature (C)')

plot_groups %+% subset( group_vals, measure == 'VWC' & season == 'spring') + ylab ( "average spring soil moisture") + xlab ( "plot group")
  
# model formula and models ----------------------------------------------------------------------- 

basic_form <- formula(v ~ (1|unique_port) + (1|date) + factor( PrecipGroup )  + Treatment*season )

m_air <- lmer(data =  df_air_temp, formula = basic_form)
m_soil <- lmer(data = df_soil_temp, formula = basic_form)

m_5cm <- lmer( data = df_5cm_soil, formula = update( basic_form, . ~ . + (1|plot) ))
m_25cm <- lmer( data = df_25cm_soil, formula = basic_form) 


summary(m_air)
summary(m_soil)
summary(m_5cm) 
summary(m_25cm)


# ------------------------------------------------------------------------------------------------- 

# make prediction df to view lmer effects 
pred_df <- df %>% select(PrecipGroup, season, Treatment) %>% distinct()

pred_df$air_temp_pred <- predict(m_air, pred_df, re.form = NA)
pred_df$VWC_5cm_pred <- predict( m_5cm, pred_df, re.form = NA )
pred_df$VWC_25cm_pred <- predict( m_25cm, pred_df, re.form = NA)

pred_stats <- pred_df %>% 
  group_by( Treatment, season ) %>% 
  summarise( air_pred = mean(air_temp_pred), VWC_5cm_pred = mean(VWC_5cm_pred), VWC_25cm_pred = mean(VWC_25cm_pred)) %>% 
  gather( type, predicted, air_pred:VWC_25cm_pred) 

pred_plot <- 
  ggplot( pred_stats, aes( x = season, y = predicted, fill = Treatment )) + 
  geom_bar( stat = 'identity', position = 'dodge') + 
  ggtitle( 'Treatment effects predicted by lmer model')

pred_stats$ylabel <- factor( pred_stats$type , labels = c('Air temperature (C)', '5 cm soil moisture', '25cm soil moisture'))

pred_plots <- pred_stats %>% group_by( type ) %>% do( p =  pred_plot %+% . + ylab(.$ylabel) ) 

pred_plots$p

# print plots in one pdf ---------------------------------------------------------------------- 

pdf( "figures/summary_of_treatment_effects.pdf" , height = 8, width = 10 ) 

print( brplt %+% subset(plot_vals, measure == 'VWC' ) + ylab( 'Soil Moisture') ) 

print( brplt %+% subset(plot_vals, measure == 'C' & depth != '25 cm deep') + ylab( 'Temperature (C)') ) 

print( brplt_tod %+% plot_T_vals + ylab( 'Temperature (C)') ) 

print ( plot_groups %+% subset( group_vals, measure == 'VWC' & season == 'spring') + ylab ( "average spring soil moisture") + xlab ( "plot group") ) 

print( pred_plots$p ) 

dev.off() 


