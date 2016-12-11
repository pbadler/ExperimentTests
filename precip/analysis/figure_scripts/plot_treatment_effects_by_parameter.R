library(stringr)
library(dplyr)
library(tidyr)
library(rstan)
library(gridExtra)
rm(list = ls())

fit_files <- dir( 'output/stan_fits', '.*climate_fit.RDS', full.names = T)
fit_files <- fit_files[- grep('best', fit_files)]

gs_fits <- fit_files [ -grep( 'recruitment', fit_files)]
r_fits <- fit_files [ grep( 'recruitment', fit_files)]

dat_files <- dir( 'data/temp_data', 'modified_.*_data_lists_for_stan.RDS', full.names = T)
treatment_stan_fit  <- dir( 'output/stan_fits', 'treatment_effects', full.names = T)
gs_treatment_fit <- treatment_stan_fit[ -grep('recruitment', treatment_stan_fit)]
r_treatment_fit <- treatment_stan_fit[ grep('recruitment', treatment_stan_fit)]

i = 2
for( i in 2:length(gs_fits)){
  
  bname <- basename(gs_fits[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
  
  dat <- readRDS( dat_files[grep( vr, dat_files)] )[[spp]]
  
  b2 <-  rstan::extract( readRDS(gs_fits[i], 'b2'))$b2
  b2 <- as.matrix( b2 ) 
  
  Chold <- as.matrix( dat$Chold )
  
  nifx <- 1:ncol(Chold)
  ifx <- grep( 'logarea', colnames(Chold))
  colnames(Chold)
  if(length(ifx) > 0 ) { 
    nifx <- nifx[-ifx]
  }
  
  Chold1 <- Chold [ , nifx, drop = F]   # climate main effects 
  Chold2 <- Chold [ , ifx , drop = F]   # climate interactions with plant size 
  
  Chold2 <- sweep( Chold2 , 2, dat$Xhold, '/') # divide by size to isolate annual climate effect 
  Chold2 <- round( Chold2 , digits = -log( .Machine$double.eps + 1, base = 10)) # round off machine error 

  alpha <- matrix( NA, nrow = nrow(b2), ncol = nrow(Chold))
  beta <- matrix( NA, nrow = nrow(b2), ncol =  nrow(Chold))
  
  for( j in 1:nrow(b2)){ 
    alpha[j, ] <- Chold1 %*% b2[ j, nifx ]
    beta[j, ]  <- Chold2 %*% b2[ j,  ifx ]
  }
  
  if ( sum(beta) != 0 ) { 
    pred_df <- data.frame( rbind( data.frame ( par = 'a (Intercept)', t(alpha) ), data.frame( par = 'b1 (size effect)', t(beta)) ) )   
  } else { 
    pred_df <- data.frame ( par = 'a (Intercept)', t(alpha) )
  }

  # make dataframe for predictions ------------------------------------
  pred_df$year <- dat$yearhold
  pred_df$Treatment <- dat$treathold
  
  rm( b2, alpha, beta, Chold, Chold1, Chold2, dat )
  
  pred_df <- unique(pred_df)
  
  pred_df <- 
    pred_df %>% 
    mutate( type = 'predicted_effect') %>% 
    gather( iteration, val , starts_with('X')) 
  
  pred_df$Treatment <- factor(pred_df$Treatment, labels = c('Control', 'Drought', 'Irrigation'))  
  
  pred_df <- 
    pred_df %>% 
    spread( Treatment, val ) %>%
    mutate( Drought = Drought - Control, Irrigation = Irrigation - Control ) 
  
  pred_df <- 
    pred_df %>% 
    gather( Treatment , val, Control, Drought, Irrigation  ) %>% 
    group_by( Treatment, par, type, iteration ) %>% 
    filter( year > 2010 ) %>%
    summarise( val = mean(val))
  
  # get observed treatment effects ---------------------------------------- # 
  bt <- rstan::extract(readRDS(gs_treatment_fit[i]), 'bt')$bt
  
  alpha <- bt[, 1:2]
  beta  <- bt[, 3:4]  
  
  obs_df <- data.frame( rbind( data.frame ( par = 'a (Intercept)', t(alpha) ), data.frame( par = 'b1 (size effect)', t(beta)) ) ) 
  
  # observation data frame   
  obs_df$Treatment <- c('Drought', 'Irrigation')
  
  obs_df <- 
    obs_df %>% 
    mutate( type = 'observed_effect') %>% 
    gather( iteration, val , starts_with('X')) 
  
  # ------------------------------------------------------------------------ # 
  pred_df <- pred_df %>% filter( Treatment != 'Control')
  
  pred_df <- as.data.frame( pred_df)
  
  df <- rbind(obs_df, pred_df )
  
  df$type <- factor( df$type)
  df$Treatment <- factor(df$Treatment ) 

  my_colors1 <- c('#1b9e77', '#d95f02', '#7570b3')
  my_colors2 <- c('black', 'orange')

  gp <- 
    ggplot(df , aes( x = val, fill = type ))  + 
    geom_density(  alpha = 0.4) + 
    geom_vline( aes(xintercept = 0), linetype = 2) +
    facet_grid( Treatment ~ . )  + 
    scale_fill_manual(values = my_colors2) +
    xlab( paste0('Treatment Effect on ', spp, ' ', vr) )
  
  p <- 
    df %>% 
    group_by( par ) %>% 
    do(p = gp %+% . + ggtitle( unique( .$par ) ) )
  
  pdf( paste( 'figures/predictions/predicted_', spp, '_', vr, '_treatment_effects.pdf' ),  width = 8, height = 6)
  
  if(length(p$p) > 1 ) { 
    print( 
      grid.arrange(p$p[[1]] + theme(legend.position = c(0.1, 0.9) ), p$p[[2]] + guides(fill = FALSE), ncol = 2, nrow = 1 )  
    )
  }else{
    print( 
      p$p[[1]]
    )
  }
    
  dev.off()  
  
  error <- 
    df %>% 
    spread(type, val ) %>% 
    mutate( prediction_error = predicted_effect - observed_effect ) 
  
  gp <- 
    ggplot( error, aes( x = prediction_error, fill = Treatment ) ) + 
    geom_density(  alpha = 0.4) + 
    geom_vline( aes(xintercept = 0), linetype = 2) +
    facet_grid( Treatment ~  . )  + 
    scale_fill_manual(values = my_colors1[ 2:3]) + 
    xlab( paste0('Predicted - observed treatment effect on ', spp, ' ', vr) )
  
  p <- 
    error %>% 
    group_by( par  ) %>% 
    do( p = gp %+% . + ggtitle(.$par ) )
  
  
  pdf( paste( 'figures/predictions/prediction_error_', spp, '_', vr, '_treatment_effects.pdf' ), width = 8, height = 6)
  
  if(length(p$p) > 1 ) { 
    print( 
      grid.arrange(p$p[[1]] + theme(legend.position = c(0.1, 0.9) ), p$p[[2]] + guides(fill = FALSE), ncol = 2, nrow = 1 )  
    )
  }else{
    print( 
      p$p[[1]]
    )
  }
  
  dev.off()  
  rm( p, gp, error, df, pred_df, obs_df)  
}
