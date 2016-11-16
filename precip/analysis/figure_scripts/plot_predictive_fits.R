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

mfiles <- dir('output/stan_fits/predictions/best_lppd', '4_predict.RDS', full.names = TRUE)

for( i in 1:length(mfiles)){ 
  
  bname <- basename(mfiles[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
  m <- mpars[3]
  lambda <- mpars[4]
  
  temp_fit <- readRDS(mfiles[i])
  
  df <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')[[spp]]
  
  # functions ---------------------------------------------------------------------------------# 
  
  plot_posterior <- function( df ) { 
    ggplot( df, aes( x = val)) + 
      geom_density() + 
      geom_vline( aes(xintercept = 0), linetype = 2, alpha = 0.7) + 
      geom_vline( aes(xintercept = lcl95), color = 'red') +
      geom_vline( aes(xintercept = ucl95), color = 'red') +
      facet_wrap( ~ var, ncol = 1 )
  }
  
  # make traceplots---------------------------------------------------------------------------------# 
  fit_ggs <- ggs(temp_fit)
  
  pdf(file.path( 'figures', 'stan_fits',  paste(spp, vr, m, 'traceplots.pdf', sep = '_')), height = 8, width = 11)  
  
  if( sum(str_detect(fit_ggs$Parameter, 'sig_a')) > 0 ){
    print( traceplot(temp_fit, c('sig_a', 'a')))
  }
  
  if( sum(str_detect(fit_ggs$Parameter, 'sig_b1')) > 0 ){
    print( traceplot(temp_fit, c('sig_b1', 'b1', 'b1_mu')))
  }

  if( sum(str_detect(fit_ggs$Parameter, 'theta')) > 0 ){
    print( traceplot(temp_fit, c('theta')))
  }
  
  if( sum(str_detect(fit_ggs$Parameter, 'tau')) > 0 ){
    print( traceplot(temp_fit, c('tau', 'tauSize')))
  }
  
  if( sum( str_detect( fit_ggs$Parameter, 'b2')) > 0  ) { 
    print( traceplot(temp_fit, c('b2')))
  }
  
  if( sum( str_detect(fit_ggs$Parameter, 'w')) > 0 ) {
    print( traceplot(temp_fit, c('w')))
  }
  
  dev.off()
  
  # plot posterior of parameter estimates for climate and crowding ----------------------------# 
  
  climate_vars <- colnames( df$C)
  
  crowding_vars <- colnames(df$W)
  
  if( sum( str_detect( fit_ggs$Parameter, 'b2')) > 0  ) { 
    b2 <- data.frame( rstan::extract( temp_fit, 'b2') ) 
    names(b2) <- climate_vars[1:ncol(b2)]
    
    b2_long <- 
      b2 %>% 
      gather( var, val, 1:ncol( b2)) %>% 
      group_by( var) %>% 
      mutate( lcl95 = quantile(val, 0.025), lcl90 = quantile(val, 0.05), med = quantile(val, 0.5), ucl90 = quantile(val, 0.95), ucl95 = quantile(val, 0.975))
    
    pdf( file.path( 'figures', 'stan_fits', paste( spp, vr, m, 'climate_effects.pdf')), height = 8, width = 11)
    
    print( plot_posterior(b2_long) + ggtitle(paste('Posterior of climate effects on', spp, vr, 'model', m)) ) 
    
    dev.off()
  } 
  
  
  if (sum(str_detect(fit_ggs$Parameter, 'w')) > 0 ) { 
    w <- data.frame(rstan::extract(temp_fit, 'w'))
    
    names(w) <- climate_vars[1:ncol(w)]
  
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
    
    pdf( file.path( 'figures', 'stan_fits', paste( spp, vr, m, 'crowding_effects.pdf')), height = 8, width = 11)
    
    print( plot_posterior(w_long) + ggtitle(paste('Posterior of crowding effects on', spp, vr, 'model', m)))
    
    dev.off()
  
  }
}
