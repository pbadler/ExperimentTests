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

mfiles <- dir('output/stan_fits', '_climate_fit.RDS', full.names = TRUE)
mfiles
for( i in 1:length(mfiles)){ 
  i = 6
  bname <- basename(mfiles[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
    
  temp_fit <- readRDS(mfiles[i])
  df <- readRDS(paste0( 'data/temp_data/modified_', vr, '_data_lists_for_stan.RDS'))[[spp]]
  
  # make traceplots---------------------------------------------------------------------------------# 

  pdf(file.path( 'figures', 'stan_fits',  paste(spp, vr, 'traceplots.pdf', sep = '_')), height = 8, width = 11)  
  
  model_pars <- temp_fit@model_pars
  
  if( 'sig_a' %in% model_pars ){
    print( traceplot(temp_fit, c('sig_a', 'a')))
  }
  
  if( 'sig_b1' %in% model_pars ){
    print( traceplot(temp_fit, c('sig_b1', 'b1', 'b1_mu')))
  }

  if( 'theta' %in% model_pars ){
    print( traceplot(temp_fit, c('theta')))
  }
  
  if( 'b2' %in% model_pars  ) { 
    print( traceplot(temp_fit, c('b2')))
  }
  
  if( 'w' %in% model_pars ) {
    print( traceplot(temp_fit, c('w')))
  }
  
  dev.off()
  
  # plot posterior of year effects ------------------------------------------------------------# 
  
  years <- unique( df$year ) 
  a <- data.frame( rstan::extract(temp_fit, 'a') ) 
  names(a)  <- years
  a$par <- 'a (intercept)'
  
  if( 'b1' %in% model_pars ){
    b1 <- rstan::extract(temp_fit, 'b1')$b1
    b1_mu <- rstan::extract(temp_fit, 'b1_mu')$b1_mu
    b1 <- sweep(b1, 1, b1_mu , '-') # subtract mean size effect to focus on deviations 
    b1 <-  as.data.frame( b1 )
    names(b1) <- years 
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
  
  pdf( file.path( 'figures', 'stan_fits', paste( spp, vr, 'year_effects.pdf', sep = '_')), height = 11, width = 8)
  
  print( plot_posterior_year_effects(year_effects_long) + ggtitle(paste('Posterior of year effects on', spp, vr, 'model')) ) 
  
  dev.off()
  
  # plot posterior of parameter estimates for climate and crowding ----------------------------# 
  
  if( is.null(colnames(df$C))){ 
    climate_vars <- names( df$Ccenter ) 
  }else{ 
    climate_vars <- colnames(df$C)
    crowding_vars <- colnames(df$W)
    if(is.null(crowding_vars)){ 
      crowding_vars <- str_extract( colnames(df$parents1), '[A-Z]{4}')
    } 
  }
   
  if( 'b2' %in% model_pars ) { 
    b2 <- data.frame( rstan::extract( temp_fit, 'b2') ) 
    names(b2) <- climate_vars[1:ncol(b2)]
        
    b2_long <- 
      b2 %>% 
      gather( var, val, 1:ncol( b2)) %>% 
      group_by( var) %>% 
      mutate( lcl95 = quantile(val, 0.025), lcl90 = quantile(val, 0.05), med = quantile(val, 0.5), ucl90 = quantile(val, 0.95), ucl95 = quantile(val, 0.975))
    
    pdf( file.path( 'figures', 'stan_fits', paste( spp, vr, 'climate_effects.pdf', sep = '_')), height = 8, width = 11)
    
    print( plot_posterior(b2_long) + ggtitle(paste('Posterior of climate effects on', spp, vr, 'model')) ) 
    
    dev.off()
  } 
  
  
  if ( 'w' %in% model_pars ) { 
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
    
    pdf( file.path( 'figures', 'stan_fits', paste( spp, vr, 'crowding_effects.pdf', sep = '_')), height = 8, width = 11)
    
    print( plot_posterior(w_long) + ggtitle(paste('Posterior of crowding effects on', spp, vr)))
    
    dev.off()
  
  }
  rm(w_long, b2_long)
}
