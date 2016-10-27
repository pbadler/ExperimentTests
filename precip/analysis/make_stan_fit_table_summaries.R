rm(list = ls() )

library(ggmcmc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(rstan)
library(gridExtra)

# input ------------------------------------------------------------------------------------# 
print(paste('Working directory: ' , getwd()))

model_table <- read.csv('output/best_WAIC_scores.csv')
model_table <- subset(model_table, vital_rate == 'recruitment')

for( do_line in 1:nrow(model_table)){ 
  
  do_model <- model_table[do_line, ]
  spp <- do_model$species
  vr <- do_model$vital_rate
  m <- do_model$model
  lambda <- do_model$lambda
  pars <- do_model$pars
  
  print(paste('print results for', spp, vr, 'model', m))
  
  temp_fit <- readRDS(file = file.path( 'output/stan_fits/predictions', paste(spp, vr, m, lambda, 4, 'predict.RDS', sep = '_')))
  
  temp_pars <-  eval ( parse ( text = as.character(pars)) ) 
  
  get_pars <- temp_pars [ !temp_pars %in% c('q_pred', 'q_pred2', 'log_lik', 'log_lik2', 'y_hat', 'y_hat2')]
  
  temp_fit_summary <- summary(temp_fit, get_pars) 
  
  pars_out <- row.names(temp_fit_summary$summary)
  
  df_summary <- data.frame(vital_rate = vr, species = spp, model = m , lambda = lambda, pars_out, temp_fit_summary$summary ) 
  
  if( do_line == 1 ) { 
    write.table(df_summary, file = 'output/model_fit_summaries.csv', append = FALSE, row.names = FALSE, sep = ',')
  }else{ 
    write.table(df_summary, file = 'output/model_fit_summaries.csv', col.names = FALSE, append = TRUE, row.names = FALSE, sep = ',')
  }
} 
