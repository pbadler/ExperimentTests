rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)

make_gg_PCA <- 
  function( spp , dl ) { 

  historical <- unique( data.frame( dl[[spp]]$C, Period = 'Historical', Treatment = 'Control', year = dl[[spp]]$year))
  modern <- unique( data.frame( dl[[spp]]$Chold, Period = 'Modern', Treatment = dl[[spp]]$treathold, year = dl[[spp]]$yearhold) ) 
  
  cvs_in <- historical[ colnames(dl[[spp]]$C)]
  cvs_out <- modern [ colnames(dl[[spp]]$Chold) ] 

  temp_pca <- princomp(cvs_in)

  dists <- apply( rbind( cvs_in, cvs_out), 1, function( x ) dist(rbind( rep(0, length(x)), x) ))

  temp_pca_both <- princomp(rbind(cvs_in, cvs_out))

  plot_df <- data.frame( temp_pca_both$scores,
            Period =c(as.character( historical$Period), 
                      as.character( modern$Period)),
            Treatment = c(historical$Treatment, modern$Treatment), 
            year = c(historical$year, modern$year), 
            dists = dists)
  
  plot_df$Treatment <- factor( plot_df$Treatment, labels = c('Control', 'Drought', 'Irrigation'))
  
  biplot(temp_pca_both)
  
  plot_df <- 
    plot_df %>% 
    mutate(center1 = mean(Comp.1[Period == 'Historical']), center2 = mean(Comp.2[Period == 'Historical']))
  
  df_arrows  <- data.frame( temp_pca_both$loadings[ , c(1:2)], x = temp_pca_both$center[1], y = temp_pca_both$center[2])
  
  names( df_arrows ) [ 1:2 ] <- c('xend', 'yend')
  
  df_arrows$xend <- df_arrows$xend*max(plot_df$Comp.1)
  df_arrows$yend <- df_arrows$yend*max(plot_df$Comp.2)
  
  gg_PCA <- 
    ggplot( plot_df, aes( x = Comp.1, y = Comp.2, color = Period )) + 
    geom_point(size = 2) + 
    geom_segment(data = df_arrows, aes (xend = xend, yend = yend, x = x, y = y), color = 1, arrow = arrow()) + 
    geom_text( data = df_arrows , aes( x = xend*1.5, y = yend*1.5, label = row.names(df_arrows)), color = 1) 
  
  return(gg_PCA)
  } 

climate <- readRDS('data/temp_data/all_clim_covs.RDS')

#lppd_scores <- read.csv('output/lppd_scores.csv')
error <- read.csv('output/climate_effect_prediction_error.csv')
datlist <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')

# base plot ------------------------------------------------------------- 

gg_PCA <- make_gg_PCA('ARTR', datlist)  
png(filename = 'figures/climate_pca_comparison.png', height = 6, width = 6, units = 'in', res = 300)
print( gg_PCA ) 
dev.off()

# plot lppd by distance -------------------------------------------------
gg_PCA <- make_gg_PCA('ARTR', datlist)

plot_df <- gg_PCA$data

MEA <- 
  error %>% 
  filter( type == 'error') %>%
  group_by(vital_rate, species, Treatment, year, type  ) %>% 
  summarise( MAE = mean(abs(val)) , ME = mean(val), sd_e = sd(val))

error_dist <- merge( plot_df, MEA , by = c('Treatment', 'year')) 

ggplot( error_dist, aes( x = dists, y = MAE)) + 
  geom_point() + 
  geom_smooth( method = 'lm', group = 1, se = F)

ggplot( error_dist, aes( x = dists, y = MAE)) + geom_point() + 
  facet_grid( ~ species ) + 
  geom_smooth( method = 'lm', group = 1, se = F)

ggplot( error_dist, aes( x = dists, y = abs( ME ) )) + geom_point() + 
  facet_grid( ~ species ) + 
  geom_smooth( method = 'lm', group = 1, se = F)

ggplot( error_dist, aes( x = dists, y = sd_e )) + geom_point() + 
  facet_grid( ~ species ) + 
  geom_smooth( method = 'lm', group = 1, se = F)



lppd_scaled_scores <- 
  lppd_scores %>% 
  filter( Period == 'Modern') %>% 
  group_by(vital_rate, species ) %>% 
  mutate( lppd_scaled = scale(lppd))

lppd_scaled_scores$year_id <- as.numeric ( factor( lppd_scaled_scores$year ) )

plot_df$treatment <- factor( plot_df$treatment, labels = c('Control', 'Drought', 'Irrigation'))

test <- subset( plot_df, Period =='Modern')


lppd_scaled_scores$lppd_scaled <- as.numeric(lppd_scaled_scores$lppd_scaled)

lppd_scaled_scores <- merge( lppd_scaled_scores, test, by = c('treatment', 'year_id'))

plot_dists <- 
  lppd_scaled_scores %>% 
  group_by(dists) %>% 
  summarise(mlppd = mean(lppd_scaled), ucl = quantile(lppd_scaled, 0.9), lcl = quantile( lppd_scaled, 0.1))

gg_dists <- 
  ggplot(lppd_scaled_scores, aes( x = dists, y = lppd_scaled  )) + 
  geom_errorbar(data = plot_dists, aes(x = dists, y = mlppd, ymin = lcl, ymax = ucl), alpha = 0.3)  + 
  geom_point( data = plot_dists, aes( x = dists, y = mlppd) )  


png(filename = 'figures/climate_dists.png', heigh = 4, width = 4, units = 'in', res = 300)
print( gg_dists + ylab( 'relative predictive score') + xlab( '"novelty"') ) 
dev.off()

