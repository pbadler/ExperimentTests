library(stringr)
library(dplyr)
library(tidyr)
library(rstan)

rm(list = ls())

fit_files <- dir( 'output/stan_fits', '*climate_fit.RDS', full.names = T)
dat_files <- dir( 'data/temp_data', 'modified_.*_data_lists_for_stan.RDS', full.names = T)
treatment_stan_fit  <- dir( 'output/stan_fits', 'treatment_effects', full.names = T)

for( i in 1:length(fit_files)){
  
  bname <- basename(fit_files[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
  
  gint_out <- rstan::extract(readRDS(fit_files[i]), 'gint_out')$gint_out
  climhat  <- rstan::extract(readRDS(fit_files[i]), 'climhat')$climhat
  
  thin <- 5
  gint_out <- gint_out[seq(1,nrow(gint_out), thin), ]
  climhat  <- climhat[seq(1,nrow(climhat), thin), ]
  
  dat <- readRDS( dat_files[grep( vr, dat_files)] )[[spp]]
  
  # make dataframe for predictions ------------------------------------
  pred_df <- data.frame(year = dat$yearhold, Group = dat$gidhold, yid = dat$yidhold, Treatment = dat$treathold, par = 'year_effect', t(climhat) ) 
  
  pred_df$logarea.t0 <- dat$Xhold
  pred_df$trackid <- dat$trackidhold
  pred_df$quad <- dat$quadhold
  
  pred_df <- 
    pred_df %>% 
    mutate( type = 'predicted_effect') %>% 
    gather( iteration, val , starts_with('X')) 
  
  # get observed year effects ---------------------------------------- # 
  yid_table <- data.frame( year = dat$yearhold, Treatment = dat$treathold, yid = dat$yidhold, Group = dat$gidhold )
  
  a <- rstan::extract(readRDS(treatment_stan_fit[i]), 'a')$a
  a     <- a[seq(1,nrow(a), thin), ]
  year_effects_a <- data.frame( yid = unique(dat$yidhold), par = 'a', t(a))
  
  treatEff <- rstan::extract(readRDS(treatment_stan_fit[i]), 'treatEff')$treatEff
  treatEff    <- treatEff[seq(1, nrow(treatEff), thin), ]
  treat_effects   <- data.frame( Group = dat$gidhold , year = dat$yearhold, Treatment = dat$treathold, yid = dat$yidhold, par = 'treatEff', t(treatEff))
  
  treat_effects$logarea.t0 <- dat$Xhold
  treat_effects$trackid <- dat$trackidhold
  treat_effects$quad   <- dat$quadhold
  
  yid_table$logarea.t0 <- dat$Xhold 
  yid_table$trackid <- dat$trackidhold
  yid_table$quad <- dat$quadhold
  
  year_effects_a <- merge( yid_table, year_effects_a, by = 'yid')
  
  if( vr != 'recruitment'){
    b1 <- rstan::extract(readRDS(treatment_stan_fit[i]), 'b1')$b1
    b1_mu <- rstan::extract(readRDS(treatment_stan_fit[i]), 'b1_mu')$b1_mu
    b1    <- b1[seq(1,nrow(b1), thin), ]
    b1_mu <- b1_mu[seq(1,nrow(b1_mu), thin) ]  
    b1 <- sweep(b1, 1, b1_mu,  '-')  
    
    year_effects_b1 <- data.frame(yid = unique(dat$yidhold), par = 'b1', t(b1))
    
    year_effects_b1 <- merge( yid_table, year_effects_b1, by = 'yid')
    
    obs_df <- rbind( year_effects_a, year_effects_b1, treat_effects)
    
    obs_df <- 
      obs_df %>% 
      mutate( type = 'observed_effect') %>% 
      gather( iteration, val, starts_with('X'))
  
    obs_df <-
      obs_df  %>% 
      spread( par, val) %>% 
      mutate( year_effect = a + b1*logarea.t0 + treatEff ) %>% 
      gather( par, val , a:year_effect)
    
  }else{ 
    
    obs_df <- rbind( unique( year_effects_a ) , unique( treat_effects )) 

    obs_df <-
      obs_df %>%
      mutate( type = 'observed_effect') %>%
      gather( iteration, val, starts_with('X'))

    obs_df <-
      obs_df %>%
      spread( par, val) %>%
      mutate( year_effect = a + treatEff ) %>%
      gather( par, val , a:year_effect)
    
  }

  # ------------------------------------------------------------------------ # 
  pred_df <- rbind(obs_df, pred_df )
  
  pred_df$Treatment <- factor(pred_df$Treatment, labels = c('Control', 'Drought', 'Irrigation'))  
  
  mean_pred_eff <-
    pred_df %>% 
    filter (par == 'year_effect') %>%
    group_by(year, Treatment, type, iteration, par ) %>% 
    summarise( val = mean( val) )

  my_colors <- c('#1b9e77', '#d95f02', '#7570b3')
  
  pdf( paste( 'figures/predictions/predicted_', spp, '_', vr, '_year_effects.pdf' ))
    
    print( 
      ggplot( subset( mean_pred_eff, Treatment == 'Control'), aes( x = val, fill = type ) ) + 
        geom_density(  alpha = 0.4) + 
        geom_vline( aes(xintercept = 0), linetype = 2) +
        facet_grid( year ~  . )  + 
        scale_fill_manual(values = my_colors) + 
        xlab( paste0('Year effects on ', spp, ' ', vr) )
    )
    
  dev.off()  
  
  mean_error <- 
    mean_pred_eff %>% 
    filter( Treatment == 'Control') %>% 
    spread(type, val ) %>% 
    mutate( prediction_error = predicted_effect - observed_effect ) 
  
  pdf( paste( 'figures/predictions/prediction_error_', spp, '_', vr, '_year_effects.pdf' ))
  
  print( 
    ggplot( mean_error, aes( x = prediction_error, fill = Treatment ) ) + 
      geom_density(  alpha = 0.4) + 
      geom_vline( aes(xintercept = 0), linetype = 2) +
      facet_grid( year ~  .  )  + 
      scale_fill_manual(values = my_colors[1]) + 
      xlab( paste0('Predicted - observed year effects on ', spp, ' ', vr) )
  )
  
  dev.off()  
  
}
