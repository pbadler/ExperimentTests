rm(list =ls())
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)

#true_dat <- read.csv('~/driversdata/data/idaho_modern/allrecords_cover.csv')
#true_dat <- read.csv('~/driversdata/data/idaho_modern/speciesData/ARTR/quadratCover.csv')

#true_dat$species <- toupper( lapply( lapply( str_split( true_dat$species, ' ' ), substr, 1,2 ) , paste0, collapse = '') ) 

g <- readRDS('output/prediction_tables/HECO_growth_3_for_cover.RDS')
s <- readRDS('output/prediction_tables/HECO_survival_2.RDS')

s <- s[ s$yid %in% g$yid, ] 

true_dat <- 
  s %>% 
  group_by( Treatment, year, quad ) %>% 
  summarise( cover = 100*(sum(exp( logarea.t0 ))/10000))

df <- data.frame( g[, c('Treatment', 'Period', 'year', 'quad', 'trackID', 'var', 'species')], logarea = as.numeric(levels(g$mean)[g$mean]), survives = as.numeric(levels(s$mean)[s$mean]))

hist( 100*(exp(df$logarea )/10000) ) 
hist( 100*(exp( s$logarea.t0)/10000))

df <- 
  df %>% 
  mutate( year = year + 1 ) %>%  # shift year up one to represent logarea.t0 to match the survival frame 
  group_by(Treatment, year, quad, var) %>% 
  summarise( cover = 100*(sum(exp(logarea)*survives))/10000)

cover <- 
  left_join(true_dat, df, by = c('Treatment', 'quad', 'year')) 

cover <- 
  cover %>% 
  gather(type, cover, cover.x, cover.y )

ggplot( subset( cover, type == 'cover.x'), aes( x = year, y = cover, color = Treatment ) ) + 
  stat_summary(fun.y = 'mean', geom = 'line') + 
  stat_summary(data = subset(cover, var == 'muhat4' & type == 'cover.y'), fun.y = 'mean', geom = 'line', linetype = 2) + 
  facet_wrap( ~ Treatment )

g_pred <- g
g_pred$year <- g$year + 1 

test <- left_join(s, g_pred , c('year', 'Treatment', 'quad', 'trackID') )

plot( test$logarea.t0.x,  as.numeric(levels(test$mean.y)[test$mean.y]))  
abline(0 , 1)
