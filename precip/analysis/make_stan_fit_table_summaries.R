rm(list = ls() )

library(ggmcmc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(rstan)
library(gridExtra)

# input ------------------------------------------------------------------------------------# 
setwd('~/Documents/ExperimentTests/precip/')
print(paste('Working directory: ' , getwd()))

model_table <- read.csv('output/best_WAIC_scores.csv')

for( do_line in 1:nrow(model_table)){ 
  
  do_model <- model_table[do_line, ]
  spp <- do_model$species
  vr <- do_model$vital_rate
  m <- do_model$model
  prior <- do_model$prior
  pars <- do_model$pars
  
  print(paste('print results for', spp, vr, 'model', m))
  
  temp_fit <- readRDS(file = file.path( 'output/stan_fits/predictions', paste(spp, vr, m, prior, 4, 'predict.RDS', sep = '_')))
  
  temp_pars <-  eval ( parse ( text = as.character(pars)) ) 
  
  get_pars <- temp_pars [ !temp_pars %in% c('mu', 'muhat', 'log_lik', 'y_hat')]
  
  temp_fit_summary <- summary(temp_fit, get_pars) 
  
  pars_out <- row.names(temp_fit_summary$summary)
  
  df_summary <- data.frame(vital_rate = vr, species = spp, model = m , prior = prior, pars_out, temp_fit_summary$summary ) 
  
  if( do_line == 1 ) { 
    write.table(df_summary, file = 'output/model_fit_summaries.csv', append = FALSE, row.names = FALSE, sep = ',')
  }else{ 
    write.table(df_summary, file = 'output/model_fit_summaries.csv', col.names = FALSE, append = TRUE, row.names = FALSE, sep = ',')
  }
} 
