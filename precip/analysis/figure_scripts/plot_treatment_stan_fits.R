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
  ggplot( df, aes( x = val)) + 
    geom_density() + 
    geom_vline( aes(xintercept = 0), linetype = 2, alpha = 0.7) + 
    geom_vline( aes(xintercept = lcl95), color = 'red') +
    geom_vline( aes(xintercept = ucl95), color = 'red') +
    facet_wrap( ~ var, ncol = 1 ) 
}

plot_posterior_year_effects <- function(df){
  ggplot( df, aes( x = val)) + 
    geom_density() + 
    geom_vline( aes(xintercept = 0), linetype = 2, alpha = 0.7) + 
    facet_grid( var ~ par)
}

# input ------------------------------------------------------------------------------------# 
setwd('~/Documents/ExperimentTests/precip/')

mfiles <- dir('output/stan_fits', '.*_treatment_fit.RDS', full.names = TRUE)
i = 7
for( i in 1:length(mfiles)){ 
  
  bname <- basename(mfiles[i])
  mpars <- unlist( str_split(bname, '_') ) 

  spp <- mpars[1]
  vr <- mpars[2]
   
  temp_fit <- readRDS(mfiles[i])
  
  df <- readRDS(paste0( 'data/temp_data/', vr, '_data_lists_for_stan.RDS'))[[spp]]
   
  model_pars <- temp_fit@model_pars
  
  # make traceplots---------------------------------------------------------------------------------# 
  
  pdf(file.path( 'figures', 'stan_fits',  paste(spp, vr, 'treatment_model_traceplots.pdf', sep = '_')), height = 8, width = 11)  
  
  if( 'sig_a' %in% model_pars ){
    print( traceplot(temp_fit, c('sig_a', 'a')))
  }
  
  if( 'sig_b1' %in% model_pars ){
    print( traceplot(temp_fit, c('sig_b1', 'b1', 'b1_mu')))
  }
  
  if( 'theta' %in% model_pars ){
    print( traceplot(temp_fit, c('theta')))
  }
  
  if( 'bt' %in% model_pars  ) { 
    print( traceplot(temp_fit, c('bt')))
  }
  
  if( 'w' %in% model_pars ) {
    print( traceplot(temp_fit, c('w')))
  }
  
  dev.off()
  
  # plot posterior of year effects ------------------------------------------------------------# 
  
  yid <- unique( data.frame( df$yidhold, df$yearhold) ) 
  
  a <- data.frame( rstan::extract(temp_fit, 'a') ) 

  names(a)  <- yid$df.yearhold
  a$par <- 'a (intercept)'
  
  if( 'b1' %in% model_pars ){
    b1_mu <- rstan::extract(temp_fit, 'b1_mu')$b1_mu
    b1 <- rstan::extract(temp_fit, 'b1')$b1
    b1 <- sweep(b1, 1, b1_mu , '-') # subtract mean size effect to focus on deviations 
    b1 <-  as.data.frame( b1 )
    names(b1) <- yid$df.yearhold 
    b1$par <- 'b1 (size effect)'
    year_effects <- rbind( a, b1)
  }else{ 
    year_effects <- a
  }
  
  year_effects_long <- 
    year_effects %>% 
    gather( var, val, - par ) %>% 
    group_by( par, var) %>% 
    mutate( lcl95 = quantile(val, 0.025), lcl90 = quantile(val, 0.05), med = quantile(val, 0.5), ucl90 = quantile(val, 0.95), ucl95 = quantile(val, 0.975))
  
  pdf( file.path( 'figures', 'stan_fits', paste( spp, vr, 'treatment_model_year_effects.pdf', sep = '_')), height = 11, width = 8)
  
  print( plot_posterior_year_effects(year_effects_long) + ggtitle(paste('Posterior of year effects on', spp, vr, 'model')) ) 
  
  dev.off()
  
  # plot posterior of parameter estimates for treatment and crowding ----------------------------# 
  
  if ( vr != 'recruitment'){ 
    treatment_vars <- c('Drought', 'Irrigation', 'Droughtxlogarea.t0', 'Irrigationxlogarea.t0')
  }else{
    treatment_vars <- c('Drought', 'Irrigation')
  }
  
  crowding_vars <- colnames(df$Whold)
  if(is.null(crowding_vars)){ 
    crowding_vars <- str_extract( colnames(df$parents1hold), '[A-Z]{4}')
  } 
  
  if( 'bt' %in% model_pars  ) { 
    bt <- data.frame( rstan::extract( temp_fit, 'bt') ) 
    names(bt) <- treatment_vars
      
    bt_long <- 
      bt %>% 
      gather( var, val, 1:ncol( bt)) %>% 
      group_by( var) %>% 
      mutate( lcl95 = quantile(val, 0.025), lcl90 = quantile(val, 0.05), med = quantile(val, 0.5), ucl90 = quantile(val, 0.95), ucl95 = quantile(val, 0.975))
    
    pdf( file.path( 'figures', 'stan_fits', paste( spp, vr, 'treatment_effects.pdf', sep = '_')), height = 8, width = 11)
    
    print( plot_posterior(bt_long) + ggtitle(paste('Posterior of treatment effects on', spp, vr, 'model')) ) 
    
    dev.off()
  } 
  
  if ('w' %in% model_pars ) { 
    w <- data.frame(rstan::extract(temp_fit, 'w'))
    rm(temp_fit)
    colnames(w) <- crowding_vars[1:ncol(w)]
    
    crowding_vars[1:ncol(w)] 
    
    if ( ncol(w) == 1 ) { 
      names(w) <- crowding_vars[ grep(crowding_vars, pattern = spp) ]
    }else{ 
      names(w) <- crowding_vars
    }
    
    w_long <- 
      w %>% 
      gather( var, val, 1:ncol(w ) ) %>% 
      group_by( var) %>% 
      mutate( lcl95 = quantile(val, 0.025), lcl90 = quantile(val, 0.05), med = quantile(val, 0.5), ucl90 = quantile(val, 0.95), ucl95 = quantile(val, 0.975))

    pdf( file.path( 'figures', 'stan_fits', paste( spp, vr, 'treatment_model_crowding_effects.pdf', sep = '_')), height = 8, width = 11)
    
    print( plot_posterior(w_long) + ggtitle(paste('Posterior of crowding effects on', spp, vr)))
    
    dev.off()
  
  }
  rm(w_long, bt_long)
}
