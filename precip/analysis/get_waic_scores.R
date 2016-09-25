###########################################################################################
#
# Read in stan output 
# calculate WAIC for each model 
# output WAIC scores in dataframe 
#
###########################################################################################

rm(list = ls())

library(loo)
library(dplyr)
library(tidyr)
library(rstan)
library(parallel)

# functions ------------------------------------------------------------------------------------


get_waic_from_file <- function( fn ) { 
  waic( 
    extract_log_lik( 
      readRDS(file = as.character( fn ))
    ) 
  )
} 

get_loo_from_file <- function( fn ) { 
  loo( 
    extract_log_lik( 
      readRDS(file = as.character( fn ))
    ) 
  )
} 


# input ------------------------------------------------------------------------------------
files <- dir('output/stan_fits/', pattern = '*.RDS', full.names = T)

file_df <- 
  data.frame(fn = files, info = gsub(pattern = '.RDS', replacement = '', basename(files)))  %>%
  separate( info , c('species', 'vital_rate', 'model', 'prior', 'chains'), sep = '_' ) %>% 
  mutate(prior = as.numeric(prior))%>%
  filter( chains > 1 ) %>% 
  arrange( species, vital_rate, model, prior) 

# run functions in parallel -----------------------------------------------------------------
file_df <- 
  file_df %>% 
  filter( species == 'POSE', model == 2 )

score_list <- mclapply(file_df$fn, get_waic_from_file, mc.cores = 3) 
waics <- unlist( lapply( score_list, function(x ) x$waic))

# score_list <- mclapply( file_df$fn, get_loo_from_file , mc.cores = 3) 
# looics <- unlist( lapply( score_list, function(x) x$looic ))
# 
# plot( unlist( lapply( score_list [ lapply( score_list, function(x) length(x)) > 1 ] , function (x ) x$looic)) )

file_df$waic <- waics 
#file_df$looic  <- looics

prior_steps <- 1:max(file_df$prior)

file_df <- merge( file_df, data.frame( prior = prior_steps, prior_stdev = seq(0.1, 1.5, length.out = max(prior_steps))))

# plot waic score vs. stdev -------------------------------------------------------------------

file_df <- 
  file_df %>% 
  group_by(species, vital_rate, model) %>% 
  mutate(min_score = ifelse(waic == min(waic), round( waic, 1 ) , ''))

score_list <- split(file_df , file_df$vital_rate)

regularization_plot <- function( x ) { 
    ggplot( subset(x, !model %in% c(1, 3)) , aes( x = prior_stdev, y = waic)) + 
      geom_point() + 
      geom_line() + 
      geom_text( aes(label = min_score), nudge_x = 0.08) + 
      facet_wrap(species ~ model, ncol = 4 ) + 
      ggtitle( paste( 'Regularization for', unique(x$vital_rate)))
}

pdf( 'figures/regularization_plots.pdf', width = 11, heigh = 8)

print( lapply( score_list, regularization_plot) ) 

dev.off()


# find best fits and rename file that has the best fit ---------------------------------------- 

best_priors <- 
  file_df %>% 
  group_by(species, vital_rate, model) %>% 
  summarise(best_prior = prior_stdev[which.min(waic)], fn = fn [ which.min(waic)], waic =  min(waic))

# output ------------------------------------------------------------------------------------

saveRDS(file_df , 'output/WAIC_scores.RDS')
saveRDS(best_priors, 'output/best_priors.RDS')
