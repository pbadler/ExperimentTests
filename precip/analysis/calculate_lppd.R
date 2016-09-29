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

# Root mean squared error ---------------------------------------------------------------------------------------------# 

compute_lppd <- function( stan_fit ) { 
  log_lik <- extract (stan_fit, "log_lik")$log_lik
  lppd <- log(colMeans(exp(log_lik)))
  lppd
} 

# input ------------------------------------------------------------------------------------# 
setwd('~/Documents/ExperimentTests/precip/')
print(paste('Working directory: ' , getwd()))

temp_fit <- readRDS(file = file.path( 'output/stan_fits/predictions/', paste(spp, vr, m, 'predictions.RDS', sep = '_')))

df <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')[[spp]]

# ---------------------------------------------------------------------------------------------------------------------

y_out <- data.frame(
  obs_id = 1:length(df$y_holdout), 
  y_holdout = df$y_holdout, 
  treatment = factor( df$treat_out, labels = c('Control', 'Drought', 'Irrigation')), 
  year = factor( df$yid_out, labels = unique(df$yid_out) + 2006 ))


# log-pointwise predictive density ------------------------------------------------------------------------------------# 

lppd <- compute_lppd(temp_fit)

y_out$lppd <- lppd  

lppd_treatment <- y_out %>% group_by(treatment) %>% summarise( mean_lppd = -mean(lppd))
lppd_year <- y_out %>% group_by(year) %>% summarise(mean_lppd = -mean(lppd))

ggplot(lppd_treatment, aes( x = treatment, y = mean_lppd)) + geom_bar(stat = 'identity')
ggplot(lppd_year, aes( x = year, y = mean_lppd )) + geom_bar(stat = 'identity')

output_lppd <- data.frame( spp = spp , vr = vr , model = m , lppd = sum(lppd))

saveRDS( output_lppd, file.path('output', paste(spp, vr, m, 'lpd_score.RDS')))
