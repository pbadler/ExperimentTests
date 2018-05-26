library(stringr)
library(dplyr)
library(tidyr)
library(rstan)
library(gridExtra)
rm(list = ls())

load('analysis/figure_scripts/my_plotting_theme.Rdata')  

fit_files <- dir( 'output/stan_fits', '.*climate_fit.RDS', full.names = T)

gs_fits <- fit_files [ -grep( 'recruitment', fit_files)]
r_fits <- fit_files [ grep( 'recruitment', fit_files)]

dat_files <- dir( 'data/temp_data', 'modified_.*_data_lists_for_stan.RDS', full.names = T)
treatment_stan_fit  <- dir( 'output/stan_fits', 'treatment_fit', full.names = T)
gs_treatment_fit <- treatment_stan_fit[ -grep('recruitment', treatment_stan_fit)]
r_treatment_fit <- treatment_stan_fit[ grep('recruitment', treatment_stan_fit)]

# first loop do all growth and survival fits ------------------------------------- # 
i = 1
for( i in 1:length(gs_fits)){
  
  bname <- basename(gs_fits[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
  
  dat <- readRDS( dat_files[grep( vr, dat_files)] )[[spp]]
  
  # extract climate parameter estimates ---------------------------------------------- # 
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
  
  Chold2 <- sweep( Chold2 , 1, dat$Xhold, '/') # divide by size to isolate annual climate effect 
  
  alpha <- matrix( NA, nrow = nrow(b2), ncol = nrow(Chold))
  beta <- matrix( NA, nrow = nrow(b2), ncol =  nrow(Chold))
  
  # calculate the overall climate effects on intercept and slope ----------- # 
  
  for( j in 1:nrow(b2)){ 
    alpha[j, ] <- Chold1 %*% b2[ j, nifx ]
    beta[j, ]  <- Chold2 %*% b2[ j,  ifx ]
  }
  
  if ( sum(beta) != 0 ) { 
    pred_df <- data.frame( rbind( data.frame ( par = 'a (Intercept)', t(alpha) ), data.frame( par = 'b1 (size effect)', t(beta)) ) )   
  } else { 
    pred_df <- data.frame ( par = 'a (Intercept)', t(alpha) )
  }

  # aggregate the climate effects by Treatment  ------------------------------------
  pred_df$year <- dat$yearhold
  pred_df$Treatment <- dat$treathold
  
  pred_df <- 
    pred_df %>% 
    mutate( type = 'predicted_effect') %>% 
    gather( iteration, val , starts_with('X'))  %>% 
    group_by( type , par, Treatment, year, iteration ) %>% 
    summarise( val = mean(val))
  
  pred_df$Treatment <- factor(pred_df$Treatment, labels = c('Control', 'Drought', 'Irrigation'))  

  pred_df <- 
    pred_df %>% 
    spread( Treatment, val ) %>%
    mutate( Drought = Drought - Control, Irrigation = Irrigation - Control ) 
  
  avg_pred_df <- 
    pred_df %>% 
    gather( Treatment , val, Control, Drought, Irrigation  ) %>% 
    group_by( Treatment, par, type, iteration ) %>% 
    filter( year > 2010 ) %>%
    summarise( val = mean(val, na.rm = T)) %>% 
    group_by( Treatment, par, type ) %>% 
    filter( !is.na(val)) %>% 
    summarise( mu = mean(val), lcl = quantile(val, 0.025), ucl = quantile(val, 0.975))

  # get observed treatment effects ---------------------------------------- # 
  bt <- read.csv(paste0( 'output/treatment_model_parameters_', spp, '_', vr, '.csv'))
  bt <- bt[ grep( '^bt', bt$X ), c(2, 5, 9)] 
  names(bt) <- c('mu', 'lcl', 'ucl')
  
  if(nrow(bt) == 4 ){ 
    obs_df <- data.frame( type = 'observed_effect', Treatment = c('Drought', 'Irrigation'), par = c('a (Intercept)', 'a (Intercept)', 'b1 (size effect)', 'b1 (size effect)'), bt )
  }else if(nrow(bt)==2){ 
    obs_df <- data.frame( type = 'observed_effect', Treatment = c('Drought', 'Irrigation'), par = c('a (Intercept)', 'a (Intercept)'), bt )
  }
  
  # bind together observed and predicted ----------------------------------- #
  
  df <- rbind( obs_df, data.frame( avg_pred_df) )
  
  df <- 
    df %>% 
    filter( Treatment != 'Control' ) %>% 
    gather( stat, val, mu:ucl) %>% 
    spread( type, val ) %>% 
    gather( type, val, observed_effect:predicted_effect ) %>% 
    spread(stat, val )
  
  # plot observed and predicted parameter estimates and confidence intervals -------------------------------- # 
  
  p1 <- 
    ggplot( df, aes( x = Treatment, y = mu, ymin = lcl, ymax = ucl, shape = type, color = Treatment)) + 
      geom_point(position = position_dodge(width = 0.3) , size = 4) + 
      geom_errorbar(position = position_dodge(width = 0.3), width = 0.3) + 
      geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.3 ) +
      ylab('Parameter estimate (+/- 95% BCI)') +
      facet_grid( par ~ . ) + 
      scale_color_manual(values = my_colors[3:4]) + 
      my_theme + 
      ggtitle(paste0( 'Treatment effects on ', spp, ' ', vr ))
    
  df$species <- spp
  df$vital_rate <- vr
  write.csv(df, paste0('output/predicted_and_observed_par_estimates_for_', spp, '_', vr, '.csv') )
    
  # plot observed and predicted effects across size ----------------------------------------------------------------- #
  bt <-  rstan::extract( readRDS(paste0( 'output/stan_fits/', spp, '_', vr, '_treatment_fit.RDS'), 'bt'))$bt
  bt <- as.matrix( bt )   
  
  tm <- dat$tmhold
  
  X <- seq(-2, 2, length.out = 100)
  
  tm <- expand.grid ( Treatment = c('Drought', 'Irrigation'), logarea.t0 = X )
  tm <- tm[ order(tm$Treatment), ]
  
  mm <- model.matrix(data = tm , ~ -1  + Treatment ) 
  
  if(ncol(bt) == 4){ 
    pred_matrix <- cbind( mm , mm*tm$logarea.t0 )
  }else{ 
    pred_matrix <- mm 
  }
  
  muhat <- matrix( ncol = nrow(pred_matrix), nrow = nrow( bt ) )
  
  for( j in 1:nrow(bt)) { i
    muhat[j , ] <- pred_matrix%*%bt[j, ]
  }
  
  mus <- colMeans(muhat)
  ci <- apply( muhat, 2, quantile, c(0.025, 0.975))
  observed_effect <- data.frame( mus, t(ci) , logarea.t0 = tm$logarea.t0, Treatment = tm$Treatment)
  names(observed_effect ) <- c('mean', 'lcl', 'ucl', 'logarea.t0', 'Treatment')

  observed_effect$type <- 'observed_effect'
  
  # get predicted effect by size -------------------------------------------------------------------------------------------# 
  
  predicted_effect <-  
    pred_df %>% 
    gather( Treatment , val, Control, Drought, Irrigation  ) %>% 
    group_by( Treatment, par, type, iteration ) %>% 
    filter( year > 2010 ) %>%
    summarise( val = mean(val, na.rm = T)) 
  
  pfx <- split(predicted_effect, predicted_effect$par )

  if(length(pfx) > 1 ) { 
    pred_bt1 <- as.matrix( pfx[[1]] %>% filter( par == 'a (Intercept)') %>% spread( Treatment, val ) %>% ungroup() %>% dplyr::select( Drought, Irrigation ))
    pred_bt2 <- as.matrix(pfx[[2]] %>% filter( par == 'b1 (size effect)') %>% spread( Treatment, val ) %>% ungroup() %>% dplyr::select( Drought, Irrigation ))
    pred_bt <- cbind(pred_bt1, pred_bt2 )
  }else if (length(pfx) == 1 ){ 
    pred_bt1 <- as.matrix( pfx[[1]] %>% filter( par == 'a (Intercept)') %>% spread( Treatment, val ) %>% ungroup() %>% dplyr::select( Drought, Irrigation ))
    pred_bt  <- pred_bt1
  } 
  
  if(ncol(bt) == 4){ 
    pred_bt <- pred_bt
  }else if(ncol(bt)==2){
    pred_bt <- pred_bt[, 1:2]
  }
  
  muhat <- matrix( ncol = nrow(pred_matrix), nrow = nrow( pred_bt ) )
  for( j in 1:nrow(pred_bt)) { 
    
    muhat[j , ] <- pred_matrix%*%pred_bt[j, ]
    
  }
  mus <- colMeans(muhat)
  ci <- apply( muhat, 2, quantile, c(0.025, 0.975))
  predicted_effect <- data.frame( mus, t(ci) , logarea.t0 = tm$logarea.t0, Treatment = tm$Treatment)
  names(predicted_effect ) <- c('mean', 'lcl', 'ucl', 'logarea.t0', 'Treatment')
  
  predicted_effect$type <- 'predicted_effect'
  
  predicted_and_observed_effects <- rbind(predicted_effect, observed_effect)
  
  p2 <- 
    ggplot (data = predicted_and_observed_effects, aes( x = logarea.t0, y = mean, 
                                                      color = Treatment, fill = Treatment, linetype = type)) + 
    geom_line() + 
    geom_hline(aes(yintercept =0 ), linetype = 2, alpha = 0.3) + 
    ylab( 'Treatment effect') + 
    xlab( 'Plant size (s.d. from mean)') + 
    scale_color_manual(values = my_colors[3:4]) + 
    my_theme


  if( nlevels( p1$data$par) == 2 ) {   
    png( paste0( 'figures/pred_v_obs_treatment_', spp, '_', vr, '.png'), width = 6, height = 4, res = 300, units = 'in')
    print(p1 )
    dev.off()
    png( paste0( 'figures/pred_v_obs_treatment_slopes_', spp, '_', vr, '.png'), width = 6, height = 4, res = 300, units = 'in')
    print(  p2 )
    dev.off( )
  }else{
    png( paste0( 'figures/pred_v_obs_treatment_', spp, '_', vr, '.png'), width = 6, height = 4, res = 300, units = 'in')
    print(p1)
    dev.off()
  }

  
  rm( p1, p2, predicted_and_observed_effects, predicted_effect, pred_matrix, observed_effect, pred_df, pfx )   
}

# end first loop,  all growth and survival fits ------------------------------------- # 


# second loop do all recruitment fits ------------------------------------------------# 

i = 1
for( i in 1:length(r_fits)){ 

  bname <- basename(r_fits[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
  
  dat <- readRDS( dat_files[grep( vr, dat_files)] )[[spp]]
  
  # extract climate parameter estimates ---------------------------------------------- # 
  b2 <-  rstan::extract( readRDS(r_fits[i], 'b2'))$b2
  b2 <- as.matrix( b2 ) 
  
  Chold <- as.matrix( dat$Chold )

  alpha <- matrix( NA, nrow = nrow(b2), ncol = nrow(Chold))

  # calculate the overall climate effects on intercept and slope ----------- # 
  for( j in 1:nrow(b2)){ 
    alpha[j, ] <- Chold %*% b2[ j,  ]
  }
  
  
  pred_df <- data.frame ( par = 'a (Intercept)', t(alpha) )
  
  # aggregate the climate effects by Treatment  ------------------------------------
  pred_df$year <- dat$yearhold
  pred_df$Treatment <- dat$treathold
  
  pred_df <- 
    pred_df %>% 
    mutate( type = 'predicted_effect') %>% 
    gather( iteration, val , starts_with('X'))  %>% 
    group_by( type , par, Treatment, year, iteration ) %>% 
    summarise( val = mean(val))
  
  pred_df$Treatment <- factor(pred_df$Treatment, labels = c('Control', 'Drought', 'Irrigation'))  
  
  pred_df <- 
    pred_df %>% 
    spread( Treatment, val ) %>%
    mutate( Drought = Drought - Control, Irrigation = Irrigation - Control ) 
  
  avg_pred_df <- 
    pred_df %>% 
    gather( Treatment , val, Control, Drought, Irrigation  ) %>% 
    group_by( Treatment, par, type, iteration ) %>% 
    filter( year > 2010 ) %>%
    summarise( val = mean(val, na.rm = T)) %>% 
    group_by( Treatment, par, type ) %>% 
    filter( !is.na(val)) %>% 
    summarise( mu = mean(val), lcl = quantile(val, 0.025), ucl = quantile(val, 0.975))
  
  # get observed treatment effects ---------------------------------------- # 
  bt <- read.csv(paste0( 'output/treatment_model_parameters_', spp, '_', vr, '.csv'))
  bt <- bt[ grep( '^bt', bt$X ), c(2, 5, 9)] 
  names(bt) <- c('mu', 'lcl', 'ucl')
  
  obs_df <- data.frame( type = 'observed_effect', Treatment = c('Drought', 'Irrigation'), par = c('a (Intercept)', 'a (Intercept)'), bt )
  
  # bind together observed and predicted ----------------------------------- #
  
  df <- rbind( obs_df, data.frame( avg_pred_df) )
  
  df <- 
    df %>% 
    filter( Treatment != 'Control' ) %>% 
    gather( stat, val, mu:ucl) %>% 
    spread( type, val ) %>% 
    mutate( predicted_effect = ifelse( is.na(predicted_effect) , 0, predicted_effect )) %>% 
    gather( type, val, observed_effect:predicted_effect ) %>% 
    spread(stat, val )
  
  # plot observed and predicted parameter estimates and confidence intervals -------------------------------- # 
  
  p1 <- 
    ggplot( df, aes( x = Treatment, y = mu, ymin = lcl, ymax = ucl, shape = type, color = Treatment)) + 
    geom_point(position = position_dodge(width = 0.3) , size = 4) + 
    geom_errorbar(position = position_dodge(width = 0.3), width = 0.3) + 
    geom_hline(aes(yintercept = 0), linetype = 2, alpha = 0.3 ) +
    ylab('Parameter estimate (+/- 95% BCI)') +
    facet_grid( par ~ . ) + 
    scale_color_manual(values = my_colors[3:4]) + 
    my_theme + 
    ggtitle(paste0( 'Treatment effects on ', spp, ' ', vr ))
  
  df$species <- spp
  df$vital_rate <- vr
  write.csv(df, paste0('output/predicted_and_observed_par_estimates_for_', spp, '_', vr, '.csv') )
  
  png( paste0( 'figures/pred_v_ob_treatment_', spp,'_', vr, '.png' ),  width = 6, height = 4, units = 'in', res = 300)
  
  print( p1 )

  dev.off()  
  
  rm( p1, predicted_effect, pred_matrix, observed_effect )   
  
}



