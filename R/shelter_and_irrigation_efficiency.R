rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(lme4)
library(zoo)

df <- readRDS('data/temp_data/decagon_data_with_rainfall_data.RDS')

# ----------------------------------------------------------------------------- 

df <- df %>% select( Treatment_label, season_label, depth_label, Treatment, PrecipGroup, plot, depth, season, unique_port, prcp_event, prerain, month, measure, datetime, date, total_rain, v, prerain, PRCP) %>% 
  filter( measure == 'VWC', month > 3, month < 11 , !is.na(prerain)) %>% 
  arrange( unique_port , datetime ) 

postrain_VWC <- df %>% 
  group_by( prcp_event) %>% 
  mutate( event_strt = min(date)) %>% 
  filter( !prerain) %>% 
  group_by(Treatment_label, season_label, depth_label, Treatment, PrecipGroup, plot, depth, season, unique_port, prcp_event, event_strt) %>% 
  summarise( cumul_rain = total_rain[which.max(v)], postrain_VWC = max(v) )

prerain_VWC <- df %>% 
  group_by( prcp_event) %>% 
  filter( prerain) %>% 
  group_by( unique_port, prcp_event) %>% 
  summarise( prerain_VWC = max(v) )

rain_effects <- left_join(postrain_VWC, prerain_VWC, by = c('unique_port', 'prcp_event'))

rain_effects <- rain_effects %>% ungroup ( ) %>% mutate( change_VWC = postrain_VWC - prerain_VWC ) 

rain_effect_plot <- ggplot( rain_effects, aes( x = cumul_rain, y = change_VWC, color = Treatment_label) ) + 
  geom_point()  + 
  geom_smooth( method = 'lm', se = FALSE) + 
  facet_wrap( ~ season_label ) + 
  ylab ( 'Change in volumetric water content (%)') + 
  xlab ( 'Cumulative rainfall during event (mm)') 

rain_effect_plots <- rain_effects %>% 
  group_by ( depth_label) %>% 
  do( p =  rain_effect_plot %+% . + ggtitle( paste0( "Effect of rain on soil moisture at ",  .$depth_label )) )

print( rain_effect_plots$p ) 

m1 <- lmer(data = subset( rain_effects, depth == '5 cm deep', season = 'spring'), change_VWC ~ prerain_VWC*cumul_rain + cumul_rain*Treatment + factor( PrecipGroup ) + (1|plot) + (1|prcp_event) )

summary(m1)

pred_df <- expand.grid( prerain_VWC = mean(rain_effects$prerain_VWC, na.rm = TRUE), 
                        cumul_rain = seq(min(rain_effects$cumul_rain, na.rm = TRUE), 
                                         max(rain_effects$cumul_rain, na.rm = TRUE), 
                                         length.out = 100), 
                        Treatment = levels( factor(rain_effects$Treatment) ) , 
                        PrecipGroup = unique( rain_effects$PrecipGroup ), 
                        depth_label = '5 cm deep', 
                        season_label = 'spring')


pred_df$change_VWC <- predict( m1, pred_df, re.form = NA)

pred_df$Treatment_label <- factor(pred_df$Treatment, levels = c('Drought', 'Control', 'Irrigation'), ordered = TRUE)


pdf( 'figures/rain_effects_plots.pdf', height = 8, width = 10 ) 

print( rain_effect_plots$p ) 

print ( rain_effect_plot %+% pred_df + ggtitle( 'Predicted rain effects at 5 cm') ) 

dev.off()


