rm(list = ls() )

library(ggmcmc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(rstan)
library(gridExtra)

# functions ---------------------------------------------------------------------------------# 

plot_posterior <- function( df ) { 
  ggplot( df, aes( x = val, color = model)) + 
    geom_density() + 
    facet_wrap(~ par , ncol = 2 ) + 
    theme_bw()
}



plot_posterior_year_effects <- function(df){
  ggplot( df, aes( x = val , color = model)) + 
    geom_density() + 
    geom_vline( aes(xintercept = 0), linetype = 2, alpha = 0.7) + 
    facet_grid( var ~ par) + 
    theme_bw()
}

# input ------------------------------------------------------------------------------------# 

mfiles <- dir('output/stan_fits', '_year_effects_fit.RDS', full.names = TRUE)
cfiles <- dir('output/stan_fits', '_climate_fit.RDS', full.names = T)
i = 1


for( i in 1:length(cfiles)){ 
  bname <- basename(cfiles[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
    
  temp_fit1 <- readRDS(cfiles[i])
  temp_fit2 <- readRDS( paste0('output/stan_fits/', spp, '_',  vr, '_year_effects_fit.RDS' ))
  
  model_pars <- temp_fit1@model_pars
  
  df <- readRDS(paste0( 'data/temp_data/', vr, '_data_lists_for_stan.RDS'))[[spp]]
  
  # plot posterior of year effects ------------------------------------------------------------# 
  
  years <- unique( df$year ) 
  a1 <- data.frame( rstan::extract(temp_fit1, 'a') ) 
  names(a1)  <- years
  a1$par <- 'a (intercept)'
  a1$model <- 'climate'
  
  a2 <- data.frame( rstan::extract(temp_fit2, 'a') ) 
  names(a2)  <- years
  a2$par <- 'a (intercept)'
  a2$model <- 'year'
  
  a <-  rbind( a1, a2)
  
  if( 'b1' %in% model_pars ){
    b1 <- rstan::extract(temp_fit1, 'b1')$b1
    b1_mu <- rstan::extract(temp_fit1, 'b1_mu')$b1_mu
    b1 <- sweep(b1, 1, b1_mu , '-') # subtract mean size effect to focus on deviations 
    b1 <-  as.data.frame( b1 )
    names(b1) <- years 
    b1$par <- 'b1 (size effect)'
    b1$model <- 'climate'
    
    b1y <- rstan::extract(temp_fit2, 'b1')$b1
    b1_muy <- rstan::extract(temp_fit2, 'b1_mu')$b1_mu
    b1y <- sweep(b1y, 1, b1_muy , '-') # subtract mean size effect to focus on deviations 
    b1y <-  as.data.frame( b1y )
    names(b1y) <- years 
    b1y$par <- 'b1 (size effect)'
    b1y$model <- 'year'

    b1 <- rbind( b1, b1y)

    year_effects <- rbind( a, b1)
    
    
  }else{ 
    year_effects <- a
  }
  
  year_effects_long <- 
    year_effects %>% 
    gather( var, val, -model,  - par ) %>% 
    group_by(model, par, var) %>% 
    mutate( lcl95 = quantile(val, 0.025), lcl90 = quantile(val, 0.05), med = quantile(val, 0.5), ucl90 = quantile(val, 0.95), ucl95 = quantile(val, 0.975))
  
  pdf( file.path( 'figures', 'stan_fits', paste( spp, vr, 'compare_year_effects.pdf', sep = '_')), height = 11, width = 8)
  
  print( plot_posterior_year_effects(year_effects_long) + ggtitle(paste('Posterior of year effects on', spp, vr, 'model')) ) 
  
  dev.off()
  
  # plot posterior of parameter estimates for climate and crowding ----------------------------# 
  
  if ( 'sig_a' %in% model_pars ) { 
    sig_a <- rstan::extract(temp_fit1, 'sig_a')$sig_a
    sig_a2 <- rstan::extract(temp_fit2, 'sig_a')$sig_a
    
    sig_a <-  data.frame( val = sig_a )
    
    sig_a$par <- 'sigma_a (intercept)'
    sig_a$model <- 'climate'
    
    sig_a2 <-  data.frame( val = sig_a2 )
    sig_a2$par <- 'sigma_a (intercept)'
    sig_a2$model <- 'year'
    
    sig_a <- rbind( sig_a, sig_a2)
    
    if( 'sig_b1' %in% model_pars ) { 
      sig_b1 <-rstan::extract(temp_fit1, 'sig_b1')$sig_b1
      sig_b12 <-rstan::extract(temp_fit2, 'sig_b1')$sig_b1
      sig_b1 <- data.frame( val  = sig_b1 )
      sig_b1$par <- 'sigma_b1 (size effect)'
      sig_b1$model <- 'climate'
    
      sig_b12 <-  data.frame( val = sig_b12 )
      sig_b12$par <- 'sigma_b1 (size effect)'
      sig_b12$model <- 'year'
    
      sig_b1 <- rbind(sig_b1, sig_b12)
      year_sd <- rbind( sig_a, sig_b1)
    }else{ 
    
      year_sd <- sig_a  
    }
    
    
  
    
    year_sd <- 
      year_sd %>% 
      group_by( model, par) %>% 
      mutate( lcl95 = quantile(val, 0.025), lcl90 = quantile(val, 0.05), med = quantile(val, 0.5), ucl90 = quantile(val, 0.95), ucl95 = quantile(val, 0.975))
    
    
    pdf( file.path( 'figures', 'stan_fits', paste( spp, vr, 'compare_year_effects_sd', sep = '_')), height = 8, width = 11)
    
    print( plot_posterior(year_sd) + ggtitle(paste('Posterior of year effects variance', spp, vr)))
    
    dev.off()
  
  }
  rm(a,year_effects, year_effects_long, year_sd)
}
