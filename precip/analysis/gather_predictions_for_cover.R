rm(list = ls())
# gather predictions rm(list = ls() )

library(dplyr)
library(tidyr)
library(stringr)
library(rstan)

# input ------------------------------------------------------------------------------------# 
setwd('~/Documents/ExperimentTests/precip/')

best_fits <- read.csv('output/WAIC_selected_models.csv')
best_fits <- subset(best_fits, vital_rate == 'growth')

best_fits$stanfile <- paste0( 'output/stan_fits/predictions/', str_replace( best_fits$fn, 'WAIC.csv', 'predict.RDS'))
best_fits$dfiles <- paste0( 'data/temp_data/', best_fits$species, '_scaled_', best_fits$vital_rate, '_dataframe.RDS')

best_fits$stanfile[1]

# get predictions -------------------------------------------------------# 
get_predictions <- function( stan_fit, y = 'muhat' ) { 
  data.frame(var = y, rstan::summary(stan_fit, y)$summary) # output summary for each data point
}
# ------------------------------------------------------------------------

for( i in 1:nrow(best_fits)){ 
  
  temp_fit <- readRDS(best_fits$stanfile[i])
  
  temp_df  <- readRDS(paste0('data/temp_data/', best_fits$species[i], '_scaled_survival_dataframe.RDS'))
  
  growth_df <- readRDS(paste0('data/temp_data/', best_fits$species[i], '_scaled_growth_dataframe.RDS'))
  
  temp_df <- temp_df[ temp_df$yid %in% growth_df$yid, ] 
  
  hold_data <- subset(temp_df, Period == 'Modern')
  
  # gather predictions ------------------------------------------------------------------------------------# 

  mu_hat1 <- get_predictions(temp_fit, 'muhat3')
  mu_hat2 <- get_predictions(temp_fit, 'muhat4')
    
  rm(temp_fit)
    
  mu_hat1 <- cbind(hold_data, mu_hat1)
  mu_hat2 <- cbind(hold_data, mu_hat2) 
    
  y_out <- rbind(mu_hat1, mu_hat2)
    
  rm(temp_df, mu_hat1, mu_hat2)

  y_out$model <- best_fits$model[i]
  y_out$vital_rate <- best_fits$vital_rate[i] 
  
  saveRDS(y_out , paste0('output/prediction_tables/', best_fits$species[i], '_', best_fits$vital_rate[i], '_', best_fits$model[i], '_for_cover.RDS'))
  
  rm(y_out)
}

