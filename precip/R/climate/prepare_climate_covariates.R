#####################################################################################
#
# Make climate variables from quarterly climate 
#
#####################################################################################

rm(list = ls() )

library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)

# ------- load files ------------------------------------------------------------------ 

quarterly_clim <- readRDS('data/temp_data/quarterly_climate.RDS')

# ------ calculate quarterly lags -----------------------------------------------------# 
#
#   Variable names follow these conventions: 
#   
#     First letter gives variable type: 
#       "P" is cumulative precipitation
#       "T" is average mean monthly temperature
#
#     Letters after the first period give the season aggregation window:
#
#       w  = winter (Q1) 
#       sp = spring (Q2)
#       su = summer (Q3)
#       f  = fall   (Q4)
#       a  = annual (Q1-4)
#
#     e.g. "P.sp" is the cumulative precipitation of the spring quarter and "P.w.sp" is
#     the cumulative precipitation of the winter and spring. 
#   
#     Number after the second period indicates the year of the transition, For example, 
#     "P.sp.1" gives the cumulative precipitation of the spring preceding the first year of the 
#     transition, whereas "T.f.w.2" gives the average temperature of the fall and winter 
#     preceding the second year of the transition. "0" refers to year before first year, 
#     i.e. "lag effect" (sensu Adler). 
#
# -------------------------------------------------------------------------------------# 
q_precip <- 
  quarterly_clim %>% 
  filter( var == 'TPCP_ttl') %>%
  group_by(Treatment) %>% 
  arrange(Treatment, year, quarter) %>%
  mutate(P.sp.1 = val, 
         P.sp.0 = lag(val, 4), 
         P.w.sp.1 = rollsum(val, 2, align = 'right', fill = NA), 
         P.w.sp.0 = lag(P.w.sp.1, 4), 
         P.su.1 = lag(val, 3), 
         P.su.0 = lag(P.su.1, 4),
         P.a.1 = rollsum(val, 4, align = 'right', fill = NA), 
         P.a.0 = lag(P.a.1, 4), 
         P.a.l = lag(P.a.0, 4)) %>% 
  filter( quarter == 'Q2') %>% # plants are measured at the end of Q2 each year 
  select( Treatment, Period, year, quarter, starts_with("P"))

q_temp <- 
  quarterly_clim %>% 
  filter( var == 'MNTM_avg') %>% 
  group_by(Treatment) %>% 
  arrange(Treatment, year, quarter) %>% 
  mutate( T.sp.1 = val, 
          T.sp.0 = lag(T.sp.1, 4),
          T.sp.l = lag(T.sp.0, 4), 
          T.w.sp.1 = rollapply(val, 2, 'mean', na.rm = TRUE, align = 'right', fill = NA), 
          T.w.sp.0 = lag(T.w.sp.1, 4),
          T.w.sp.l = lag(T.w.sp.0, 4),
          T.su.1 = lag(val, 3),
          T.su.0 = lag(T.su.1, 4)) %>% 
  filter( quarter == 'Q2') %>% 
  select( Treatment, Period, year, quarter, starts_with("T"))

allClim <- 
  q_precip %>% 
  left_join ( q_temp, by = c('Treatment', 'Period', 'quarter', 'year')) %>% 
  arrange( Treatment, year) 

allClim$year <- allClim$year - 1 # adjust to match Peter's assignment of year 0 as the reference year in demographic data sets

# -- calculate interactions --------------------------------------------------------------#

# make interactions: 
# test$TxP.sp.1 <- test$P.sp.1*test$T.sp.1
# test$TxP.sp.2 <- test$P.sp.2*test$T.sp.2


# ---- output ----------------------------------------------------------------------------# 

saveRDS( data.frame( allClim ) , 'data/temp_data/all_clim_covs.RDS')
