rm(list =ls())
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(boot)
library(rstan)
# functions ---------------------------------------------------------------------------- 

#

setwd('~/Documents/ExperimentTests/precip/')

my_colors <- c('#1b9e77', '#d95f02', '#7570b3')

years <- expand.grid(year = 1925:2017, Treatment = c(1:3), stat = c('true_cov', 'pred_cover'))
years$Period[ years$year > 2006 ] <- 'Modern'
years$Period[ years$year <= 1960 ] <- 'Historical'

# output lists 
survives <- list(NA)
size     <- list(NA)
recruits <- list(NA)

species_list <- c('ARTR', 'HECO', 'POSE', 'PSSP')

ylims <- list( c(0,40), c(0,7.5), c(0,7.5), c(0,7.5))
i = 1

for( i in 1:4) {  
  spp   <- species_list[i]  
  
  sdl   <- readRDS('data/temp_data/modified_survival_data_lists_for_stan.RDS')[[i]]
  rdl   <- readRDS('data/temp_data/modified_recruitment_data_lists_for_stan.RDS')[[i]]
  gdl   <- readRDS('data/temp_data/modified_growth_data_lists_for_stan.RDS')[[i]]
  
  m <- dir('output/stan_fits', paste0( spp, '.*_climate_fit.RDS'), full.names = TRUE)
  m1 <- dir('output/stan_fits', paste0( spp, '.*year_effects_fit.RDS'), full.names = TRUE)
  
  smu <- summary(readRDS(m[3]), c('muhat'))$summary[, 1]
  plot( sdl$Xhold , smu ) 
  
  fit <- readRDS(m[1])
  fity <- readRDS(m1[1])
  
  gmu <- summary(fit, c('mu'))$summary[,1]
  df1 <- data.frame(obs = gdl$obs_id, Y = gdl$Y, mu = gmu, type = '1')
  plot(df1$mu, df1$Y) 
  abline( 0, 1, col = 'red')
  
  gmu <- summary(fit, c('muhat'))$summary[,1]
  df2 <- data.frame(obs = gdl$obs_idhold, Y = gdl$Yhold, mu = gmu, 'type' = 'hold')
  plot( df1$mu, df1$Y , col = 'gray')
  points(df2$mu, df2$Y) 
  abline( 0, 1, col = 'red')
  
  df <- rbind( df1, df2)
  
  gmu <- summary(fit, c('muhat3'))$summary[,1]
  df3 <- data.frame(obs = gdl$obs_id3, mu3 = gmu)
  
  df <- df %>% 
    spread(type, mu) 
  
  df <- merge( df, df3)
  
  plot(df$mu3, df$Y) 
  abline( 0, 1, col = 'red')
  points(df$hold, df$Y, col = 'blue')
  
  gmu <- summary( fity, 'muhat')$summary[,1]
  df2y <- data.frame( obs = gdl$obs_idhold, Y = gdl$Yhold, mu = gmu, type  = 'hold')
  plot( df2y$mu, df2y$Y)
  points(df$hold, df$Y, col = 'gray')
  abline(0,1)
  
  test <- read.csv('data/temp_data/ARTR_growth_and_survival_cleaned_dataframe.csv')
  
  plot(test$logarea.t0, test$logarea.t1, col = 'gray')
  points(data = subset(test, Period == 'Modern'), logarea.t1 ~ logarea.t0)
  m1 <- lm(data = subset(test, Period == 'Historical'), logarea.t1 ~ logarea.t0)
  abline(m1)
  m2 <- lm(data = subset(test, Period == 'Modern'), logarea.t1 ~ logarea.t0)
  abline(m2, col = 'blue')
  abline(0,1, lty = 2)
  
  
  
  gmu <- summary(readRDS(m[1]), c('muhat'))$summary[,1]
  
  pdf( paste0( 'figures/predictions/', spp  , '_predicted_cover.pdf' ), height = 8, width = 8) 
  
  print( 
    ggplot( subset( plot_cover[[i]], Period == "Historical"), aes( x = year, y =  val, color = Treatment, linetype = stat)) +
      geom_line() +
      scale_color_manual(values = my_colors ) + 
      ylim( ylims[[i]]) + 
      scale_x_continuous()
  )
  
  print( 
    ggplot( subset( plot_cover[[i]], Period == "Modern"), aes( x = year, y =  val, color = Treatment, linetype = stat)) +
      geom_line() +
      scale_color_manual(values = my_colors ) + 
      ylim( ylims[[i]] ) + 
      scale_x_continuous()
  )
  
  dev.off()
} 

saveRDS(size , 'output/predicted_size.RDS')
saveRDS(survives, 'output/predicted_survival.RDS')
