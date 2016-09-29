rm(list = ls() )

library(ggmcmc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(rstan)
library(gridExtra)

# arguments --------------------------------------------------------------------------------# 

args <- commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 3){ 
  stop('####### Incorrect number of arguments supplied ####### \n
       ####### Arguments required:  
       #######  spp: species name ("ARTR" "HECO" "PSSP" or "POSE") 
       #######  vr:  vital rate ("growth" "survival" or "recruitment")
       #######  m:  model number: (1, 2, 3, 4, 5)')
}else if (length(args) == 3){
  
  # ---Set working directory, species, vital rate, model number, and number of chains -----------------------------#
  args <- commandArgs(trailingOnly = TRUE)
  
  spp <- args[1]
  vr <- args[2]
  m <- args[3]
  
  print(paste('plot results for', spp, vr, 'model', m))
  
}
# # # for testing
# spp <- 'POSE'
# vr <- 'growth'
# m <- 2


compute_dev <- function( stan_fit, y_obs ) { 
  
  y_hat <- extract (stan_fit, "y_hat")$y_hat
  dev <- colMeans(y_hat) - y_obs 
  
  dev
  
} 


# input ------------------------------------------------------------------------------------# 
setwd('~/Documents/ExperimentTests/precip/')
print(paste('Working directory: ' , getwd()))

temp_fit <- readRDS(file = file.path( 'output/stan_fits/predictions', paste(spp, vr, m, 'predictions.RDS', sep = '_')))

df <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')[[spp]]

# plot observed - predicted  ------------------------------------------------------------

y_out <- data.frame(
  obs_id = 1:length(df$y_holdout), 
  y_holdout = df$y_holdout, 
  treatment = factor( df$treat_out, labels = c('Control', 'Drought', 'Irrigation')), 
  year = factor( df$yid_out, labels = unique(df$yid_out) + 2006 ))

muhat <- extract( temp_fit, 'muhat')$muhat 
y_hat <- extract( temp_fit, 'y_hat')$y_hat

mu_diff <- colMeans( sweep(muhat, MARGIN = 2, y_out$y_holdout, FUN = '-'))
y_hat_diff <- colMeans( sweep( y_hat, MARGIN = 2, y_out$y_holdout, FUN = '-'))

y_out$mu_diff <- mu_diff 
y_out$y_hat_diff <- y_hat_diff 

g1 <- 
  ggplot(y_out, aes(x = mu_diff)) + 
  geom_histogram() + 
  geom_vline( aes(xintercept = 0), linetype = 2, alpha = 0.5) + 
  ggtitle( paste( 'Predicted - observed for', spp, vr, 'model', m, '\n(mu_hat - observed)'))

g2 <- 
  ggplot(y_out, aes(x = y_hat_diff)) + 
  geom_histogram() + 
  geom_vline( aes(xintercept = 0), linetype = 2, alpha = 0.5) + 
  ggtitle( paste( 'Predicted - observed for', spp, vr, 'model', m, '\n(y_hat - observed)'))

pdf( file.path( 'figures', paste(spp, vr, m, 'all_predicted_minus_observed.pdf')),  width = 11, height = 8 )

print( grid.arrange(g1, g2, nrow = 2) ) 

dev.off()

# 
g_treatments <- 
  ggplot( y_out, aes( x = mu_diff)) + 
  geom_density() + 
  facet_wrap(~ treatment, ncol = 1 ) +  
  geom_vline( aes(xintercept = 0), linetype = 2, alpha = 0.5) + 
  ggtitle( paste( 'Predicted - observed for', spp, vr, 'model', m, '\n(mu_hat - observed)'))

g_years <- 
  ggplot( y_out, aes( x = mu_diff)) + 
  geom_density() + 
  facet_wrap( ~ year, ncol = 1 ) + 
  geom_vline( aes( xintercept = 0 ) , linetype = 2, alpha = 0.5 ) + 
  ggtitle( paste( 'Predicted - observed for', spp, vr, 'model', m, '\n(mu_hat - observed)'))


pdf( file.path( 'figures/predictions', paste(spp, vr, m, 'predicted_minus_obs_treatments.pdf')), width = 8, height = 11)
print( g_treatments)
dev.off()

pdf( file.path( 'figures/predictions', paste(spp, vr, m, 'predicted_minus_obs_years.pdf')), width = 8, height = 11)
print( g_years)
dev.off()

# Root mean squared error ---------------------------------------------------------------------------------------------# 

y_out$deviations <- compute_dev(temp_fit, y_out$y_holdout)

RMSE_by_treatment <- y_out %>% group_by( treatment) %>% summarise( RMSE = sqrt( sum(deviations^2)/n() ))
RMSE_by_year <- y_out %>% group_by( year, treatment ) %>% summarise( RMSE = sqrt(sum(deviations^2)/n()))

ggplot( RMSE_by_year, aes( x = year, y = RMSE)) + geom_bar(stat = 'identity') + facet_wrap( ~ treatment)
ggplot( RMSE_by_treatment, aes( x = treatment, y = RMSE )) + geom_bar(stat = 'identity')




