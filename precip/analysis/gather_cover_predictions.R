rm(list = ls())
# gather predictions rm(list = ls() )

library(stringr)
library(rstan)

# input ------------------------------------------------------------------------------------# 
model_table <- read.csv('output/best_WAIC_scores.csv')

mfiles <- dir('output/stan_fits/predictions', '.*_growth_.*_4_predict.RDS', full.names = TRUE)
mfiles <- head(mfiles, 2)

# get predictions -------------------------------------------------------# 

get_predictions <- function( stan_fit, y = 'muhat' ) { 
  cbind(var = y, rstan::summary(stan_fit, y)$summary) # output summary for each data point
}

# ---------------------------------------------------------------------------------------------------------------------

for( i in 1:length(mfiles)){ 
  
  bname <- basename(mfiles[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
  m <- mpars[3]
  lambda <- mpars[4]
  
  temp_fit <- readRDS(mfiles[i])
  
  temp_df  <- readRDS(paste0('data/temp_data/', spp, '_scaled_survival_dataframe.RDS'))
  
  growth_df <- readRDS(paste0('data/temp_data/', spp, '_scaled_growth_dataframe.RDS'))
  
  temp_df <- temp_df[ temp_df$yid %in% growth_df$yid, ] 
  
  hold_data <- subset(temp_df, Period == 'Modern')
  
  # gather predictions ------------------------------------------------------------------------------------# 
  
  mu_hat3 <- get_predictions(temp_fit, 'muhat3')
  mu_hat4 <- get_predictions(temp_fit, 'muhat4')
    
  mu_hat3 <- cbind(hold_data, mu_hat3)
  mu_hat4 <- cbind(hold_data, mu_hat4) 
    
  y_out <- rbind(mu_hat3, mu_hat4)
    
  rm(temp_fit, temp_df, growth_df, mu_hat3, mu_hat4)
  
  y_out$model <- m
  y_out$vital_rate <- vr 
  
  saveRDS(y_out , paste0('output/prediction_tables/', spp, '_', vr, '_', m, '_for_cover.RDS'))
  rm(y_out)
}