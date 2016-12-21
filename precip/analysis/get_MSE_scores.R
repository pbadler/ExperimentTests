rm(list =ls())
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(boot)
library(rstan)
# functions ---------------------------------------------------------------------------- 

#
load('analysis/figure_scripts/my_plotting_theme.Rdata')

setwd('~/Documents/ExperimentTests/precip/')

years <- expand.grid(year = 1925:2017, Treatment = c(1:3), stat = c('true_cov', 'pred_cover'))
years$Period[ years$year > 2010 ] <- 'Modern'
years$Period[ years$year <= 2010 ] <- 'Historical'

# output lists 
survives <- list(NA)
size     <- list(NA)
recruits <- list(NA)

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')

i = 1
out <- out2 <- list()

for( i in 1:4) {  
  spp   <- species_list[i]  
  
  sdl   <- readRDS('data/temp_data/modified_survival_data_lists_for_stan.RDS')[[i]]
  rdl   <- readRDS('data/temp_data/modified_recruitment_data_lists_for_stan.RDS')[[i]]
  gdl   <- readRDS('data/temp_data/modified_growth_data_lists_for_stan.RDS')[[i]]
  
  m <- dir('output/stan_fits', paste0( spp, '.*_climate_fit.RDS'), full.names = TRUE)
  m1 <- dir('output/stan_fits', paste0( spp, '.*year_effects_fit.RDS'), full.names = TRUE)
  
  sdf <- data.frame(sdl[ c('treathold', 'Yhold', 'gidhold', 'yearhold')] )
  rdf <- data.frame(rdl[ c('treathold', 'Yhold', 'gidhold', 'yearhold')] )
  gdf <- data.frame(gdl[ c('treathold', 'Yhold', 'gidhold', 'yearhold')] )
  
  sdf$climate <- summary(readRDS(m[3]), c('muhat'))$summary[, 1]
  rdf$climate <- summary(readRDS(m[2]), c('lambda_pred'))$summary[, 1]
  gdf$climate <- summary(readRDS(m[1]), c('muhat'))$summary[, 1]
  
  sdf$year_effects <- summary(readRDS(m1[3]), c('muhat'))$summary[, 1]
  rdf$year_effects <- summary(readRDS(m1[2]), c('lambda_pred'))$summary[, 1]
  gdf$year_effects <- summary(readRDS(m1[1]), c('muhat'))$summary[, 1]
  
  sdf$vital_rate <- 'survival'
  gdf$vital_rate <- 'growth'
  rdf$vital_rate <- 'recruitment'
  
  df <- rbind( sdf, gdf, rdf )
  
  names(df) <- c('Treatment', 'Yhold', 'Group', 'year', 'climate', 'year_effects', 'vital_rate')
  df$Treatment <- factor(df$Treatment, labels = c('Control', 'Drought', 'Irrigation'))
  
  df <- 
    df %>% 
    gather( model, mu , climate:year_effects )
  
  df$SE <- (df$Yhold - df$mu)^2  
  df$species <- spp 
  
  overall_MSE <- 
    df %>% 
    group_by( species, vital_rate, model ) %>% 
    summarise( MSE = mean(SE))
  
  treatment_MSE <-   
    df %>% 
    group_by( species, vital_rate, Treatment, model ) %>% 
    summarise( MSE = mean(SE))
  
  out[[i]] <- overall_MSE
  out2[[i]] <- treatment_MSE
} 

out <- do.call( rbind, out )
out2 <- do.call(rbind, out2)

write.csv(out, 'output/all_MSE_scores.csv', row.names = F)
write.csv(out2, 'output/treatment_MSE_scores.csv', row.names = F)
