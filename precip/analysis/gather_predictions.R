rm(list = ls())
# gather predictions rm(list = ls() )

library(dplyr)
library(tidyr)
library(stringr)
library(rstan)

# input ------------------------------------------------------------------------------------# 
setwd('~/Documents/ExperimentTests/precip/')

model_table <- read.csv('output/best_WAIC_scores.csv')

best_fits <- model_table %>% group_by(species, vital_rate ) %>% filter( waic == min(waic))

best_fits$stanfile <- paste0( 'output/stan_fits/predictions/', str_replace( best_fits$fn, 'WAIC.csv', 'predict.RDS'))
best_fits$dfiles <- paste0( 'data/temp_data/', best_fits$species, '_scaled_', best_fits$vital_rate, '_dataframe.RDS')

mfiles <- dir('output/stan_fits/predictions', '4_predict.RDS', full.names = TRUE)

# get predictions -------------------------------------------------------# 
get_predictions <- function( stan_fit, y = 'muhat' ) { 
  data.frame(var = y, rstan::summary(stan_fit, y)$summary) # output summary for each data point
}
# ------------------------------------------------------------------------

for( i in 1:nrow(best_fits)){ 
  
  temp_fit <- readRDS(best_fits$stanfile[i])
  
  temp_df  <- readRDS(best_fits$dfiles[i])
  
  hold_data <- subset(temp_df, Period == 'Modern')
  
  # gather predictions ------------------------------------------------------------------------------------# 

  if(best_fits$vital_rate[i] == 'survival'){ 
    
    mu_hat1 <- get_predictions(temp_fit, 'muhat')
    mu_hat2 <- get_predictions(temp_fit, 'muhat2')
    
    rm(temp_fit)

    mu_hat1 <- cbind(hold_data, mu_hat1)
    mu_hat2 <- cbind(hold_data, mu_hat2) 
    
    y_out  <- rbind(mu_hat1, mu_hat2)
    
    rm(temp_df, mu_hat1, mu_hat2)
  }else if(best_fits$vital_rate[i] == 'growth'){ 
    
    mu_hat1 <- get_predictions(temp_fit, 'muhat')
    mu_hat2 <- get_predictions(temp_fit, 'muhat2')
    
    rm(temp_fit)
    
    mu_hat1 <- cbind(hold_data, mu_hat1)
    mu_hat2 <- cbind(hold_data, mu_hat2) 
    
    y_out <- rbind(mu_hat1, mu_hat2)
    
    rm(temp_df, mu_hat1, mu_hat2)
  }else if(best_fits$vital_rate[i] == 'recruitment'){ 
    
    lambda_hat1 <- get_predictions(temp_fit, 'lambda_pred')
    lambda_hat2 <- get_predictions(temp_fit, 'lambda_pred2')
    
    rm(temp_fit)
    
    lambda_hat1 <- cbind(hold_data, lambda_hat1)
    lambda_hat2 <- cbind(hold_data, lambda_hat2) 
    
    y_out <- rbind(lambda_hat1, lambda_hat2)
    
    rm(temp_df, lambda_hat1, lambda_hat2)
  }
  
  y_out$model <- best_fits$model[i]
  y_out$vital_rate <- best_fits$vital_rate[i] 
  
  saveRDS(y_out , paste0('output/prediction_tables/', best_fits$species[i], '_', best_fits$vital_rate[i], '_', best_fits$model[i], '.RDS'))
  
  rm(y_out)
}

