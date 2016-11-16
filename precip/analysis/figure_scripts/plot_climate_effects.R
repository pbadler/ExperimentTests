library(stringr)
library(dplyr)
library(tidyr)
library(rstan)

rm(list = ls())

fit_files <- dir( 'output/stan_fits/predictions/', '4_predict.RDS', full.names = T)
dat_files <- dir( 'data/temp_data/', 'data_lists_for_stan.RDS', full.names = T)
treatment_stan_fit  <- dir( 'output/stan_fits', 'treatment_effects', full.names = T)
for( i in 1:length(fit_files)){
  
  bname <- basename(fit_files[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]
  m <- mpars[3]
  lambda <- mpars[4]
  
  gint_out <- rstan::extract(readRDS(fit_files[i]), 'gint_out')$gint_out
  climhat  <- rstan::extract(readRDS(fit_files[i]), 'climhat')$climhat
  
  thin <- 10
  gint_out <- gint_out[seq(1,nrow(gint_out), thin), ]
  climhat  <- climhat[seq(1,nrow(climhat), thin), ]
  
  dat <- readRDS( dat_files[grep( vr, dat_files)] )[[spp]]
  
  # make dataframe for predictions ------------------------------------
  pred_df <- data.frame(year = dat$yearhold, Group = dat$gidhold, yid = dat$yidhold, Treatment = dat$treathold, t(climhat) ) 
  
  pred_df <- 
    pred_df %>% 
    mutate( type = 'predicted_effect') %>% 
    gather( iteration, val , starts_with('X')) 
  
  # get observed treatment effects ---------------------------------------- # 
  gint <- rstan::extract(readRDS(treatment_stan_fit[i]), 'gint')$gint
  treatEff <- rstan::extract(readRDS(treatment_stan_fit[i]), 'treatEff')$treatEff
  gint <- gint[seq(1,nrow(gint), thin), ]
  treatEff  <- treatEff[seq(1,nrow(treatEff), thin), ]
  
  obs_df <- data.frame(year = dat$yearhold, Group = dat$gidhold, yid = dat$yidhold, Treatment = dat$treathold, t(treatEff) ) 

  obs_df <- 
    obs_df %>% 
    mutate( type = 'observed_effect') %>% 
    gather( iteration, val, starts_with('X'))
  
  # ------------------------------------------------------------------------ # 
  
  pred_df <- rbind(obs_df, pred_df )
  pred_df$Treatment <- factor(pred_df$Treatment, labels = c('Control', 'Drought', 'Irrigation'))  
  
  mean_pred_eff <-
    pred_df %>% 
    filter (year > 2010 & Group == 1) %>%
    group_by(Treatment, type, iteration ) %>% 
    summarise( val = mean( val) ) 
  
  mean_pred_eff <- 
    mean_pred_eff %>% 
    spread(Treatment, val ) %>% 
    group_by( type ) %>%
    mutate( mControl = mean(Control)) %>%
    mutate( Drought  = Drought - mControl, Irrigation = Irrigation - mControl) %>%
    gather( Treatment, val , Control:Irrigation ) %>% 
    group_by( type, Treatment ) %>% 
    mutate( mean_effect = mean(val))

  my_colors <- c('#1b9e77', '#d95f02', '#7570b3')
  
  pdf( paste( 'figures/predictions/predicted_', spp, '_', vr, '_treatment_effects.pdf' ))
    
    print( 
      ggplot( mean_pred_eff %>% filter(Treatment != 'Control'), aes( x = val, fill = type ) ) + 
        geom_density(  alpha = 0.4) + 
        geom_vline( data = mean_pred_eff %>% filter( Treatment != 'Control' , type == 'observed_effect'), aes(xintercept = mean_effect), linetype = 2, alpha = 0.8) +
        #geom_vline( aes(xintercept = 0), linetype = 2) +
        facet_grid( Treatment ~  .  )  + 
        scale_fill_manual(values = my_colors) + 
        xlab( 'Treatment Effect on Growth')
    )
    
  dev.off()  
}
