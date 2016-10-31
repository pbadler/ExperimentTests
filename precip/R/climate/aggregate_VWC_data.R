#######################################################################################
#
# Setup seasonal SoilWAT variables for demographic rate models 
#
#######################################################################################

rm(list = ls()) 

library( ggplot2 ) 
library(tidyr)
library(dplyr)
library(lme4)
library(zoo)
library(stringr)

# make time periods --------------------------------------------------------------------

p1 <- data.frame( Period = 'Modern', year = 2007:2016)
p2 <- data.frame( Period = 'not monitored', year = 1958:2006)
p3 <- data.frame( Period = 'Historical', year = 1925:1957)
periods <- data.frame( rbind( p1, p2, p3 )) 

# ----- read in data --------------------------------------------------------------------#

VWC <- readRDS('data/temp_data/daily_swVWC_treatments.RDS')

# ---process dates----------------------------------------------------------------------#

any(is.na(VWC$month))
VWC$date
VWC$season
VWC$year <- as.numeric( VWC$year  ) 

# set-up aggregate seasonal variables for model ----------------------------------------#

df <- 
  VWC %>% 
  ungroup() %>% 
  mutate( water_year = year + lag_year ) %>% 
  mutate( quarter = cut(month, 4, labels = paste0('Q', 1:4))) %>%
  dplyr::select(year, quarter, month, year, season, season_label, precip_seasons, water_year, Treatment, layer, date, VWC, VWC_raw)

df <- df %>% left_join(periods, by = 'year')

df <- 
  df %>% 
  filter( !( Treatment != 'Control' & year < 2012) ) # remove Drought and Irrigation prior to experiment 

# ---------- annual soil moisture -------------------------------------------------#
annual_VWC <- 
  df %>% 
  group_by( Period, year, layer, Treatment ) %>%
  summarise (avg_VWC = mean(VWC_raw, na.rm = TRUE))

# ---------- seasonal soil moisture  -----------------------------------------------#
seasonal_VWC <- 
  df %>% 
  mutate(year = ifelse(month == 12 , year + 1, year  )) %>% # account for December 
  group_by(Period, year, season_label, Treatment, layer) %>% 
  summarise( l0 = mean(VWC_raw, na.rm = TRUE) )

# ------------ aggregate VWC with Treatment effects by quarter ---------------#
df[is.na(df$quarter), ]

quarterly_VWC <-
  df %>% 
  group_by( Period, Treatment, layer, year, quarter ) %>%                       
  summarise( avg = mean(VWC_raw, na.rm = TRUE)) 


# -------- output -----------------------------------------------------------------------------#

saveRDS( seasonal_VWC, 'data/temp_data/seasonal_VWC.RDS')
saveRDS( quarterly_VWC, 'data/temp_data/quarterly_VWC.RDS')
saveRDS( annual_VWC, 'data/temp_data/annual_VWC.RDS')

