# 
#  plot WAIC scores by prior 
# 
# 
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(zoo)

loess_predict <- 
  function( waic , log_lambda, error = FALSE ) { 
    test_fit <- loess(formula = waic ~ log_lambda)
    test_predict <- predict(test_fit, se = TRUE)
    test_predict$fit
    test_predict$se.fit
    if (error ) { 
      return( test_predict$se.fit )
    }else if(!error){ 
      return( test_predict$fit)
    }
    
  }

all_waics <- read.csv('output/WAIC_scores.csv')

load('data/temp_data/master_list.Rdata')
model_table <- read.csv('data/temp_data/model_table.csv')

# regularization based on Gerber et al. 2015 ---------------------------------------------------------------------# 
nlambda <- master_list$nlambda
lambda.set <- exp(seq(-5, 15, length=nlambda))
sd_vec <- sqrt(1/lambda.set) # use sd for stan normal distribution 
# ----------------------------------------------------------------------------------------------------------------# 

all_waics <- 
  all_waics %>% 
  mutate( clean_fn = str_replace(fn, pattern = '_WAIC\\.csv$', replacement = '')) %>%
  separate(clean_fn , into = c('species', 'vital_rate', 'model', 'prior', 'chains'), sep = '_') %>% 
  mutate( model = as.numeric(model), prior = as.numeric( prior))

waics <- 
  merge(all_waics, model_table, by = c('model', 'species', 'vital_rate', 'prior')) %>% 
  mutate( log_lambda = log(lambda.set[as.numeric(prior)]) ) %>% 
  mutate( prior_sd = sd_vec[as.numeric(prior)]) %>% 
  group_by(species, vital_rate, model) %>% 
  arrange( vital_rate, species, model , as.numeric(prior)) %>% 
  filter( n() > 1 ) %>% 
  mutate( mean_WAIC = mean(waic), 
          sd_WAIC = sd(waic), 
          outlier = ifelse(abs(waic - mean_WAIC) > 400 , TRUE, FALSE ), 
          min_waic = min(waic[!outlier])) %>% 
  group_by( species, vital_rate, model, outlier) %>% 
  mutate( lfit = ifelse( outlier, NA, loess_predict(waic , log_lambda )), 
          lfit.se = ifelse( outlier, NA, loess_predict(waic, log_lambda, TRUE)), 
          waic_diff = waic - min_waic )

gg_waic_base <- 
  ggplot( waics, aes( x = log_lambda, y = waic_diff )) + 
  geom_text( aes( x = log_lambda, y = waic_diff, label = prior))  

gg_waic <- 
  gg_waic_base +  
  geom_line(aes(x = log_lambda, y = lfit - min_waic, group = 1))

gg_no_outlier <- 
  waics %>% 
  filter (!outlier ) %>% 
  group_by(species, vital_rate, model) %>% 
  filter( n() > 1 ) %>% 
  do(gg = gg_waic_base %+% . + 
       geom_smooth(aes( group = 1), color = 'blue',  se = FALSE) + 
       ggtitle(paste('WAIC by regularization plot for', unique(.$species), unique(.$vital_rate), 'model', unique(.$model)) ))


gg <- waics %>% 
  group_by(species, vital_rate, model) %>% 
  filter( n() > 1 ) %>% 
  do(gg = gg_waic %+% . + 
       ggtitle(paste('WAIC by regularization plot for', unique(.$species), unique(.$vital_rate), 'model', unique(.$model)) ))


# output ----------------------------------------------------------------------------------------

pdf('figures/plot_WAIC_by_lambda_no_outliers.pdf', height = 8 , width = 8 )

print( gg_no_outlier$gg ) 

dev.off()


pdf('figures/plot_WAIC_by_lambda.pdf', height = 8 , width = 8 )

print( gg$gg ) 

dev.off()

# save outliers to be re-run ----------------------------------------------------------------------# 

outliers <- 
  waics %>% 
  filter( outlier ) 

print( outliers ) 

write.csv(outliers, 'output/outlier_runs.csv', row.names = FALSE)

# save lowest WAIC models 

best_fits <- 
  waics %>% 
  group_by( vital_rate, species, model ) %>% 
  filter( lfit == min(lfit[!is.na(lfit)] )) 

no_reg_fits <- 
  all_waics %>% 
  filter( model %in% c(1, 3)) 

best_fits <- bind_rows ( best_fits %>% select_(.dots = names( all_waics)), no_reg_fits )
best_fits <- left_join( best_fits, model_table)

write.csv(best_fits , 'output/best_WAIC_scores.csv')