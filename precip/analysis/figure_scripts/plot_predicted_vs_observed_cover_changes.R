# 
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
rm(list = ls() ) 

load('analysis/figure_scripts/my_plotting_theme.Rdata')

species_list <- c('ARTR', 'HECO', 'POSE' , 'PSSP')

cover_threshholds <- c(2, 0, 0.25, 0.25)
x_adj <- c(0.75, 0.75, 0.45, 0.45)
grd <- data.frame(species = species_list, cover_threshholds = cover_threshholds, x_adj = x_adj , model = c('climate'))
grd2 <- grd
grd2$model <- 'year_effects'
grd <- rbind(grd,grd2)
i = 1
out <- labels <- stats <-  list(NA)

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
  
  df$Period <- ifelse(df$year > 2011, 'Modern', 'Historical')
  
  df <- 
    df %>% 
    filter( obs*100/10000 > cover_threshholds, 
            pred*100/10000 > cover_threshholds )  # only use plots with greater than cover threshold 
  
  allcombos <- expand.grid(Treatment = c('Control', 'Drought', 'Irrigation'), year = c(1925:2016))
  
  df <- subset( df , Period == 'Modern' ) # confine it to the experimental years only 
  
  overall_MSE <- 
    df %>% ungroup %>% summarise( MSE = mean((observed-predicted)^2))
  
  overall_cor <- 
    df %>% ungroup %>% summarise( cor = cor(predicted, observed))
  
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
    scale_color_manual(values = my_colors[2:4]) + 
    scale_alpha_manual(values = c(1, 1, 1)) + 
    my_theme +
    ggtitle(paste(spp))
    # coord_fixed()
  
  pdf(paste0( 'figures/predictions/', spp, '_', model , '_model_predicted_and_observed_population_growth_rates.pdf'), height = 8, width = 8)
  print( pts )
  dev.off()
  
  overall_stats <- data.frame(species = spp, model = paste0( model, ' model') )
  overall_stats$MSE <- as.numeric(overall_MSE)
  overall_stats$cor <- as.numeric(overall_cor)
  
  label_df$species <- spp
  label_df$model <- paste0( model, ' model' )
  
  df$species <- spp
  df$model <- paste0( model, ' model')
  
  labels[[i]] <- label_df
  out[[i]] <- df 
  stats[[i]] <- overall_stats
  
  rm(df)
  
}
 
tail(out)

out <- do.call(rbind, out)
label_df <- do.call(rbind, labels )
stats <- do.call(rbind, stats)


i = 1
for( i in 1:length(species_list)) { 
  
  spp <- species_list[i]
  
  png(paste0('figures/', spp, '_predicted_pgr_comparison.png'), width = 8, height = 10, res = 300, units = 'in')
  print( 
  ggplot( subset( out, species == spp), aes(x = predicted, y = observed, color = Treatment )) +
    geom_point() +
    geom_smooth(method = 'lm', se = F, alpha = 0.2, linetype = 2, color = 1, size = 1) +
    geom_text( data = subset(label_df, species == spp), aes( x = pos.x, y = pos.y , color = NULL , label = unique(label)), hjust = 1, vjust = -0.5 , show.legend = F) +
    geom_abline(aes(intercept = 0, slope = 1), color = 'gray') + 
    facet_grid(  Treatment ~ model  )+
    xlab( 'Annual population growth rate predicted') +
    ylab( 'Annual population growth rate observed') +
    # scale_y_continuous(limits = c(-ylim, ylim)) +
    # scale_x_continuous(limits = c(-ylim, ylim)) +
    scale_color_manual(values = my_colors[2:4]) +
    my_theme +
    ggtitle(species_names[[i]])
  )
  dev.off()
}

# make stats table 
stats$model <- str_replace(stats$model, '_', ' ')
stats <- stats %>% gather( stat , val, MSE:cor)


fit_table <- 
  stats %>% 
  spread(model, val) %>% 
  mutate(diff = `climate model` - `year effects model`) %>% 
  mutate( improved = ifelse(diff > 0 & stat == 'cor', '***', '')) %>% 
  mutate( improved = ifelse(diff < 0 & stat == 'MSE', '***', improved))

fit_table <- fit_table %>% dplyr::select(species, stat, `year effects model` , `climate model` , diff , improved )

corxt <- xtable(fit_table, caption = 'MSE of predicted log cover changes and correlations between log cover changes predicted and observed. Predictions for the cover changes in the experimental plots were generated either from the year effects or the climate models. Instances where the climate model made better predictions than the year effects model are indicated with the "***". ARTR = \\textit{A. tripartita}, HECO = \\textit{H. comata}, POSE = \\textit{P. secunda}, PSSP = \\textit{P. spicata}.',
       label = 'table:corPGR')

print(corxt, 'manuscript/pgr_predictions.tex', type = 'latex', caption.placement ="top", table.placement = 'H')

# ----------- all treatments and species ---------------
out$model <- str_replace(out$model, '_', ' ')

out$model <- factor(out$model, labels = c('climate model', 'no climate model'))

out <- out[complete.cases(out),]


levels(out$model)

scores_by_treatment <- 
  out %>% 
  group_by( Treatment, species, model)  %>% 
  summarise( cor = cor( predicted , observed ), MSE = mean((predicted-observed)^2)) %>% 
  gather( stat, val, cor:MSE) %>% 
  spread(model, val) %>% 
  ungroup() %>% 
  mutate(diff = `climate model` - `no climate model`) %>% 
  mutate( improved = ifelse( stat == 'cor' & diff > 0 , '***', '')) %>% 
  mutate( improved = ifelse( stat == 'MSE' & diff < 0 , '***', improved))

fit_table2 <- 
  scores_by_treatment %>% dplyr::select(species, Treatment, stat, `no climate model` , `climate model` , diff , improved ) %>% arrange( species, Treatment, stat )

corxt2 <- xtable(fit_table2, caption = 'MSE of predicted log cover changes and correlations between log cover changes predicted and observed. Predictions for the cover changes in the experimental plots were generated either from the year effects or the climate models. Instances where the climate model made better predictions than the year effects model are indicated with the "***". ARTR = \\textit{A. tripartita}, HECO = \\textit{H. comata}, POSE = \\textit{P. secunda}, PSSP = \\textit{P. spicata}.',
                label = 'table:corPGR')

print(corxt2, 'manuscript/pgr_predictions_by_treatment.tex', type = 'latex', caption.placement ="top", table.placement = 'H')

