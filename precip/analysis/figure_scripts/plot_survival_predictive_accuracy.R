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

mfiles <- dir('output/stan_fits/predictions', '.*survival.*4_predict.RDS', full.names = TRUE)
dffiles <- dir('data/temp_data', '.*survival.*cleaned_dataframe.RDS', full.names = T)
for( i in 1:length(mfiles)){ 
  
  bname <- basename(mfiles[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
  m <- mpars[3]
  lambda <- mpars[4]
  
  temp_fit <- readRDS(mfiles[i])
  df <- readRDS(dffiles[i])
  
  thin <- 10
  mu <- rstan::extract(temp_fit, 'mu')$mu
  mu <-  mu[ seq(1, nrow(mu), thin),  ]
  muhat <- rstan::extract(temp_fit, 'muhat')$muhat
  muhat <-  muhat[ seq(1, nrow(muhat), thin),  ]
  # sigma <- rstan::extract(temp_fit, 'sigma')$sigma
  # sigma <-  sigma[ seq(1, nrow(sigma), thin),  ]
  # sigmahat <- rstan::extract(temp_fit, 'sigmahat')$sigmahat
  # sigmahat <-  sigmahat[ seq(1, nrow(sigmahat), thin),  ]
  # 
  rm(temp_fit)
  
  mu <- cbind(mu, muhat)
  # sigma <- cbind(sigma, sigmahat)
  # pred <- matrix( rnorm(length(mu), mu, sigma), dim(mu))
  
  diff <- sweep( t(mu), 1, df$Y, '-')
  
  rm(mu)
  
  diff_long <- 
    data.frame(df, diff) %>% 
    gather( observation, val , X1:X100) 
  
  rm(df, diff)
  
  diff_long_period <- 
    diff_long %>% 
    group_by( Period ) %>% 
    mutate( ucl95 = quantile(val,  0.975), lcl95 = quantile(val, 0.025))

  diff_long_treatment <- 
    diff_long %>% 
    group_by( Period, Treatment ) %>% 
    mutate( ucl95 = quantile(val,  0.975), lcl95 = quantile(val, 0.025))
  
  diff_long_year <- 
    diff_long %>% 
    group_by( year ) %>% 
    mutate( ucl95 = quantile(val,  0.975), lcl95 = quantile(val, 0.025))
  
  diff_long_year_treatment <- 
    diff_long %>% 
    group_by( year , Treatment  ) %>% 
    mutate( ucl95 = quantile(val,  0.975), lcl95 = quantile(val, 0.025))
  
  bp_hist <-  
    ggplot( diff_long_period , aes( x = val)) + 
      geom_histogram() + 
      geom_vline(aes( xintercept = 0), lty = 2 , alpha = 0.7) + 
      geom_vline(aes( xintercept = ucl95), color = 'red')+ 
      geom_vline(aes( xintercept = lcl95), color = 'red') + 
      xlab ( 'Surival probability predicted - survival outcome')
    
  # bp_dens <-  
  #   ggplot( diff_long_period , aes( x = val)) + 
  #   geom_density() + 
  #   geom_vline(aes( xintercept = 0), lty = 2 , alpha = 0.7) + 
  #   geom_vline(aes( xintercept = ucl95), color = 'red')+ 
  #   geom_vline(aes( xintercept = lcl95), color = 'red') 
  # 
  #   
  pdf(file.path( 'figures', 'stan_fits',  paste(spp, 'survival_predictive_accuracy.pdf', sep = '_')), height = 8, width = 11)  
  

    print(  bp_hist + facet_grid( Period ~ . ))
    print(  bp_hist %+% subset(diff_long_year_treatment, Period == 'Modern') + facet_grid( year ~ Treatment ))
    print(  bp_hist %+% subset(diff_long_year, Period == 'Modern') + facet_grid( year ~ . ) )
    print(  bp_hist %+% subset(diff_long_treatment, Period == 'Modern') + facet_grid(Treatment ~  .  ))
    
    rm(bp_hist)
    
    # print(  bp_dens + facet_grid( Period ~ . ))
    # print(  bp_dens %+% subset(diff_long_year_treatment, Period == 'Modern') + facet_grid( year ~ Treatment ))
    # print(  bp_dens %+% subset(diff_long_year, Period == 'Modern') + facet_grid( year ~ . ) )
    # print(  bp_dens %+% subset(diff_long_treatment, Period == 'Modern') + facet_grid(Treatment ~  .  ))
    # 
  dev.off()
  
}
