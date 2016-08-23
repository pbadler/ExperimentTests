# make season and time of data dataframes 

rm(list = ls()) 
library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(lme4)
library(zoo)


season <- data.frame ( month = 1:12, season = c('winter', 'winter', 'spring', 'spring', 'spring', 'summer', 'summer', 'summer', 'fall', 'fall', 'fall', 'winter'))

season$season_label <- factor( season$season, levels = c('spring', 'summer', 'fall', 'winter'), order = TRUE)

season$precip_seasons <- factor(c(rep('cool', 5), rep('warm', 5), rep('cool', 2)))
season$lag_year <- 0 
season$lag_year[ season$month > 10 ] <- 1 

tod <- data.frame( hour = 1:24, tod = cut(1:24, breaks = c(0, 6, 19, 24)) )
levels(tod$tod ) <- c('night', 'day', 'night')

saveRDS(season, 'data/temp_data/season.RDS')

saveRDS(tod, 'data/temp_data/tod.RDS')

