# 
library(dplyr)
library(tidyr)
library(stringr)

rm(list = ls() ) 

load('analysis/figure_scripts/my_plotting_theme.Rdata')

species_list <- c('ARTR', 'HECO', 'POSE' , 'PSSP')

cover_threshholds <- c(2, 0, 0.5, 1)
x_adj <- c(0.75, 0.75, 0.45, 0.45)
grd <- data.frame(species = species_list, cover_threshholds = cover_threshholds, x_adj = x_adj , model = c('climate'))
grd2 <- grd
grd2$model <- 'year_effects'
grd <- rbind(grd,grd2)


for(i in 1:nrow(grd)){ 
  
  spp <- grd$species[i]
  model <- grd$model[i]
  cover_threshholds <- grd$cover_threshholds[i]
  x_adj <- grd$x_adj[i]
  # get last year of cover which is not in the survival dataframe ------------------------------------- #   
  
  oc <- read.csv(paste0('output/ibm/simulations/', spp, '_observed_cover_per_quadrat.csv'))
  
  m1 <- dir( 'output/ibm/simulations', as.character(spp) , full.names = T)
  pc <- readRDS(  paste0('output/ibm/simulations/', as.character(spp), '_one_step_ahead_', as.character(model),'_model_quadrat_cover.RDS'))
  
  # ------------------------------------------------------------------------------------------------------ #
  oc$obs <- oc$area
  pc$pred <- pc$area
  df <- merge(oc[, c('year', 'Treatment', 'quad', 'obs')], pc[, c('year', 'Treatment', 'simulation', 'quad', 'pred')], by = c('year', 'Treatment', 'quad' ))
  
  pred_pgr <- 
    df %>% 
    mutate( pgr = log(pred) - log(obs)) %>%
    filter( is.finite(pgr)) %>% 
    group_by(year, Treatment, quad) %>% 
    summarise( predicted = mean(pgr), lcl = quantile(pgr, 0.25), ucl = quantile(pgr, 0.75), pred = mean(pred))
  
  pred_pgr$year = pred_pgr$year + 1 # predicted pgr for next year 
  
  # get observed pgr --------------------------------------------------------------------- # 
  
  all_years <- expand.grid( year = min(oc$year):max(oc$year), quad = unique(oc$quad))
  all_years <- merge( all_years, unique(oc[, c('quad', 'Treatment')] ) )
  
  oc <- merge(all_years, oc, by = c('year', 'quad', 'Treatment'), all.x = T )
  
  obs_pgr <- 
    oc %>% 
    group_by( Treatment , quad ) %>% 
    arrange(Treatment, quad,  year ) %>% 
    mutate( obs_lag = lag( obs, 1 )) %>% 
    mutate( observed = log(obs)- log(obs_lag))  %>% 
    filter( is.finite(observed) ) %>% 
    dplyr::select(year, quad, Treatment, obs, observed ) 
  
  # get mean squared error ---- 
  df <- merge(obs_pgr, pred_pgr)
  
  df$Period <- ifelse(df$year > 2000, 'Modern', 'Historical')
  T2 <- as.character( df$Treatment)
  T2 <- ifelse(df$Period == 'Historical' , 'Historical', T2)
  df$Treatment <- T2
  df$Treatment <- factor( df$Treatment, levels = c('Historical', 'Control', 'Drought', 'Irrigation'), ordered = T)
  
  
  df <- 
    df %>% 
    filter( obs*100/10000 > cover_threshholds, pred*100/10000 > cover_threshholds )  # only use plots with greater than cover threshold 
  
  allcombos <- expand.grid(Treatment = c('Historical', 'Control', 'Drought', 'Irrigation'), year = c(1925:2016))
  
  MSE <- 
    df %>% 
    group_by(Treatment ) %>% 
    summarise( MSE = mean( (observed - predicted)^2 ) )
  
  cor <- 
    df %>% 
    group_by( Treatment ) %>% 
    summarise( cor = cor(predicted, observed))
  
  xlim <- max( df$predicted )
  ylim <- min(  df$observed)
  
  
  label_df <- merge(cor, MSE)
  
  label_df <- 
    label_df %>% 
    mutate( pos.x = x_adj*xlim, pos.y = -Inf ) %>% 
    #mutate( pos.x = , pos.y = 0.8*ylim ) %>% 
    mutate(cor = round(cor, 2), MSE = round(MSE, 2)) %>% 
    mutate( label = paste0('r=', cor, '\n', 'MSE=', MSE))
  
  df <- merge(allcombos, df, all.x = T)
  
  pts <- 
    ggplot( df, aes( x = predicted, y = observed, color = Treatment) ) + 
    geom_point(aes(alpha = Treatment)) +
    geom_smooth(method = 'lm', se = F, alpha = 1, linetype = 1, color = 1, size = 1) + 
    geom_text(data = label_df, aes( x = pos.x, y = pos.y , color = NULL , label = label), hjust = 1, vjust = -0.5 , show.legend = F) +   
    facet_wrap( ~ Treatment  )+ 
    xlab( 'Annual population growth rate predicted') + 
    ylab( 'Annual population growth rate observed') + 
    # scale_y_continuous(limits = c(-ylim, ylim)) + 
    # scale_x_continuous(limits = c(-ylim, ylim)) + 
    scale_color_manual(values = my_colors) + 
    scale_alpha_manual(values = c(0.4, 1, 1, 1)) + 
    my_theme +
    ggtitle(paste(spp))
    # coord_fixed()
  
  pdf(paste0( 'figures/predictions/', spp, '_', model , '_model_predicted_and_observed_population_growth_rates.pdf'), height = 8, width = 8)
  print( pts )
  dev.off()
}

