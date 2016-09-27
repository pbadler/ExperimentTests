# 
#  plot WAIC scores by prior 
# 
# 
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

waics <- readRDS('output/WAIC_scores.RDS')

# regularization based on Gerber et al. 2015 ---------------------------------------------------------------------# 
nlambda <- 30
lambda.set <- exp(seq(-5, 15, length=nlambda))
sd_vec <- sqrt(1/lambda.set) # use sd for stan normal distribution 
# ----------------------------------------------------------------------------------------------------------------# 

waics <- 
  waics %>% 
  mutate( clean_fn = str_replace(fn, pattern = '^\\./', replacement = '')) %>% 
  mutate( clean_fn = str_replace(clean_fn, pattern = '\\.RDS$', replacement = '')) %>%
  separate(clean_fn , into = c('species', 'vital_rate', 'model', 'prior', 'chains'), sep = '_') %>% 
  mutate( prior_sd = sd_vec[as.numeric(prior)])

gg_waic <- ggplot( waics, aes( x = prior_sd, y = waic ) ) + geom_point() + geom_smooth() 

gg <- waics %>% group_by(species, vital_rate, model) %>% filter( n() > 1 ) %>% do(gg = gg_waic %+% . + ggtitle(paste('WAIC by sd plot for', unique(.$species), unique(.$vital_rate), 'model', unique(.$model)) ))

pdf('figures/plot_WAIC_by_sd.pdf', height = 8 , width = 8 )

print( gg$gg ) 

dev.off()
