library(stringr)
library(dplyr)
library(tidyr)
library(rstan)

rm(list = ls())

fit_files <- dir( 'output/stan_fits', '*climate_fit.RDS', full.names = T)
dat_files <- dir( 'data/temp_data', 'modified_.*_data_lists_for_stan.RDS', full.names = T)
treatment_stan_fit  <- dir( 'output/stan_fits', 'treatment_effects', full.names = T)

thin <- 5

for( i in 1:length(fit_files)){
  
  bname <- basename(fit_files[i])
  mpars <- unlist( str_split(bname, '_') ) 
  
  spp <- mpars[1]
  vr <- mpars[2]

  b2 <-  rstan::extract( readRDS(fit_files[i], 'b2'))$b2
  
  if(length(dim(b2)) == 1) {
    b2  <- b2[ seq( 1, length(b2), thin) ]
  }else{
    b2 <- b2[seq(1, nrow(b2), thin), ]
  }

  dat <- readRDS( dat_files[grep( vr, dat_files)] )[[spp]]
  Chold <- dat$Chold
  
  if ( vr == 'recruitment' ) { 
    Chold1 <- Chold   # climate main effects 
    alpha  <- matrix( NA, nrow = nrow(b2), ncol = nrow( Chold1 ))
    
    for( j in 1:nrow(b2)){ 
      alpha[j, ] <- Chold1 %*% b2[ j,  ]
    }
    
    pred_df <- data.frame(  par = 'a (Intercept)', t(alpha) )
    
  }else{ 
    
    if ( length(dim(Chold)) == 2 ){ 
      
      ifx <- grep( 'logarea', colnames(Chold))
      nifx <- (1:ncol(Chold))[-ifx]
      
      Chold1 <- Chold [ , nifx ]   # climate main effects 
      Chold2 <- Chold [ , ifx ]    # climate interactions with plant size 
      
      alpha <- matrix( NA, nrow = nrow(b2), ncol = nrow(Chold))
      beta <- matrix( NA, nrow = nrow(b2), ncol =  nrow(Chold))
      
      for( j in 1:nrow(b2)){ 
        
        if( !is.null( dim(Chold1) ) ){ 
          alpha[j, ] <- Chold1 %*% b2[ j, nifx ]
        }else{ 
          alpha[j, ] <- Chold1 * b2[j, nifx]
        }
        if( !is.null( dim(Chold2) ) ){ 
          beta[j, ]  <- ( Chold2 / dat$Xhold ) %*% b2[ j,  ifx ]
        }else{ 
          beta[j, ]  <- ( Chold2 / dat$Xhold ) * b2[ j,  ifx ] 
        }
      }
       
      pred_df <- data.frame( rbind( data.frame ( par = 'a (Intercept)', t(alpha) ), data.frame( par = 'b1 (size effect)', t(beta)) ) )   
      
    }else{ 
      
      alpha <- Chold%*%t(b2) 
      pred_df <- data.frame ( par = 'a (Intercept)', alpha  )
    }
  }
  
  # make dataframe for predictions ------------------------------------
  pred_df$yid <- dat$yidhold
  pred_df$Treatment <- dat$treathold
  pred_df$quad <- dat$quadhold
  pred_df$year <- dat$yearhold
  pred_df$Group <- dat$gidhold
  
  pred_df <- 
    pred_df %>% 
    mutate( type = 'predicted_effect') %>% 
    gather( iteration, val , starts_with('X')) 
  
  # get observed treatment effects ---------------------------------------- # 
  
  gint <- rstan::extract(readRDS(treatment_stan_fit[i]), 'gint')$gint
  bt <- rstan::extract(readRDS(treatment_stan_fit[i]), 'bt')$bt
  treatEff <- rstan::extract(readRDS(treatment_stan_fit[i]), 'treatEff')$treatEff
  gint <- gint[seq(1,nrow(gint), thin), ]
  treatEff  <- treatEff[seq(1,nrow(treatEff), thin), ]
  bt        <- bt[seq(1,nrow(bt), thin), ]
  
  tm <- dat$tmhold [ , -1] 
  
  if( vr != 'recruitment') { 
    alpha <- matrix( NA, nrow = nrow(bt), ncol = nrow(tm))
    beta  <- matrix( NA, nrow = nrow(bt), ncol = nrow(tm))
    
    for( j in 1:nrow(bt)){ 
      alpha[j, ] <- tm %*% bt[ j, 1:2] 
      beta[ j, ] <- tm %*% bt[ j, 3:4]
    }
    obs_df <- data.frame( rbind( data.frame ( par = 'a (Intercept)', t(alpha) ), data.frame( par = 'b1 (size effect)', t(beta)) ) )   
    
  }else{
    
    alpha <- matrix( NA, nrow = nrow(bt), ncol = nrow(tm))
    
    for( j in 1:nrow(bt)){ 
      alpha[j, ] <- tm %*% bt[ j, 1:2] 
    }
    
    obs_df <-  data.frame ( par = 'a (Intercept)', t(alpha) )
    
  }
  
  obs_df$yid <- dat$yidhold
  obs_df$Treatment <- dat$treathold
  obs_df$quad <- dat$quadhold
  obs_df$year <- dat$yearhold
  obs_df$Group <- dat$gidhold
  
  obs_df <- 
    obs_df %>% 
    mutate( type = 'observed_effect') %>% 
    gather( iteration, val , starts_with('X')) 
  
  # ------------------------------------------------------------------------ # 
  df <- rbind(obs_df, pred_df )
  
  df$Treatment <- factor(df$Treatment, labels = c('Control', 'Drought', 'Irrigation'))  
  
  mean_eff <-
    df %>% 
    filter (year > 2010 & Group == 1) %>%
    group_by(Treatment, type, iteration, par  ) %>% 
    summarise( val = mean( val) )

  mean_eff <- 
    mean_eff %>% 
    spread(Treatment, val ) %>% 
    mutate( Drought  = Drought - Control, Irrigation = Irrigation - Control) %>%
    gather( Treatment, val , Control:Irrigation )

  my_colors <- c('#1b9e77', '#d95f02', '#7570b3')
  
  pdf( paste( 'figures/predictions/predicted_', spp, '_', vr, '_treatment_effects.pdf' ))
    
    print( 
      ggplot( mean_eff %>% filter(Treatment != 'Control'), aes( x = val, fill = type ) ) + 
        geom_density(  alpha = 0.4) + 
        geom_vline( aes(xintercept = 0), linetype = 2) +
        facet_grid( Treatment ~  par  )  + 
        scale_fill_manual(values = my_colors) + 
        xlab( paste0('Treatment Effect on ', spp, ' ', vr) )
    )
    
  dev.off()  
  
  mean_error <- 
    mean_eff %>% 
    filter( Treatment != 'Control') %>% 
    spread(type, val ) %>% 
    mutate( prediction_error = predicted_effect - observed_effect ) 
  
  pdf( paste( 'figures/predictions/prediction_error_', spp, '_', vr, '_treatment_effects.pdf' ))
  
  print( 
    ggplot( mean_error, aes( x = prediction_error, fill = Treatment ) ) + 
      geom_density(  alpha = 0.4) + 
      geom_vline( aes(xintercept = 0), linetype = 2) +
      facet_grid( Treatment ~  par  )  + 
      scale_fill_manual(values = my_colors[2:3]) + 
      xlab( paste0('Predicted - observed treatment effect on ', spp, ' ', vr) )
  )
  
  dev.off()  
  
}
