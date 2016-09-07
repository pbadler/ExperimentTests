#####################################################################################
#
# Merge growth data with climate data; prepare for STAN
#
#####################################################################################

rm(list = ls() )

library(dplyr)
library(tidyr)
library(ggplot2)

# ------- load files ------------------------------------------------------------------ 
growth_data <- dir('data/temp_data/', pattern = 'growth.RDS', full.names = TRUE)

spp_names <- regmatches(growth_data, m = gregexpr( pattern = '([A-Z]{4})', growth_data )) 

clim_m <- readRDS('data/temp_data/monthly_climate.RDS')
seasonal_clim <- readRDS('data/temp_data/seasonal_climate.RDS')
quarterly_clim <- readRDS('data/temp_data/quarterly_climate.RDS')

# ------ calculate quarterly lags -----------------------------------------------------# 
#
#   Variable names follow these conventions: 
#   
#     First letter gives variable type: 
#       "P" is cumulative precipitation
#       "T" is average mean monthly temperature
#
#     Number after the first period gives aggregation window. e.g. "P.4" represents the 
#     cumulative precipitation of the previous four quarters  
#   
#     Number after the second period indicates the either the first year "1" or second 
#     year "2" of the transition  of the transition. For example, "P.4.1" gives the 
#     cumulative precipitation of the four quarters preceding the first year of the 
#     transition, whereas "T.8.2" gives the average temperature of the 8 quarters 
#     preceding the second year of the transition. 
#
# -------------------------------------------------------------------------------------# 

q_precip <- 
  quarterly_clim %>% 
  filter( var == 'TPCP_ttl') %>% 
  arrange(year, quarter) %>%
  mutate(P.1.2 = val, 
         P.1.1 = lag(val, 4), 
         P.2.2 = rollsum(val, 2, align = 'right', fill = NA), 
         P.2.1 = lag(P.2.2, 4), 
         P.4.2 = rollsum(val, 4, align = 'right', fill = NA), 
         P.4.1 = lag(P.4.2, 4), 
         P.8.2 = rollsum(val, 8, align = 'right', fill = NA), 
         P.8.1 = lag(P.8.2, 4)) %>% 
  filter( quarter == 'Q2') %>% 
  select( period, year, quarter, treatment, starts_with("P"))


q_temp <- 
  quarterly_clim %>% 
  filter( var == 'MNTM_avg') %>% 
  arrange( year, quarter) %>% 
  mutate( T.1.2 = val, 
          T.1.1 = lag(T.1.2, 4), 
          T.2.2 = rollapply(val, 2, 'mean', na.rm = TRUE, align = 'right', fill = NA), 
          T.2.1 = lag(T.2.2, 4), 
          T.4.2 = rollapply(val, 4, 'mean', na.rm = TRUE, align = 'right', fill = NA), 
          T.4.1 = lag(T.4.2, 4), 
          T.8.2 = rollapply(val, 8, 'mean', na.rm = TRUE, align = 'right', fill = NA), 
          T.8.1 = lag(T.8.2, 4)) %>% 
  filter( quarter == 'Q2') %>% 
  select( period, year, quarter, treatment, starts_with("T"))


allClim <- 
  q_precip %>% 
  left_join ( q_temp, by = c('period', 'year', 'quarter', 'treatment')) %>% 
  arrange( year, quarter, treatment )

# -- merge with all growth records -------------------------------------------------------#

allD <- lapply( growth_data, readRDS)

allD <- lapply(allD, function(x, ... ) { merge ( x, ...  ) }, y = allClim, by = 'year')  
                

# ---- output ----------------------------------------------------------------------------# 

saveRDS( allD, 'data/temp_data/all_growth_combined.RDS')

saveRDS(allClim, 'data/temp_data/all_clim_vars.RDS')



