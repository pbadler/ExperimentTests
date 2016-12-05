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
  function( lppd , lambda, error = FALSE ) { 
    test_fit <- loess(formula = lppd ~ lambda)
    test_predict <- predict(test_fit, se = TRUE)
    test_predict$fit
    test_predict$se.fit
    if (error ) { 
      return( test_predict$se.fit )
    }else if(!error){ 
      return( test_predict$fit)
    }
    
  }


all_lppds <- read.csv('output/WAIC_scores_compiled.csv')

all_lppds <- all_lppds %>% filter( type == 'out_of_sample') %>% dplyr::select(-type)

load('data/temp_data/master_list.Rdata')
model_table <- read.csv('data/temp_data/model_table.csv')

# ----------------------------------------------------------------------------------------------------------------# 

all_lppds <- 
  all_lppds %>% 
  mutate( clean_fn = str_replace(fn, pattern = '_WAIC\\.csv$', replacement = '')) %>%
  separate(clean_fn , into = c('species', 'vital_rate', 'model', 'lambda', 'chains'), sep = '_') %>% 
  mutate( model = as.numeric(model), lambda = as.numeric( lambda))

lppds <- 
  merge(all_lppds, model_table, by = c('model', 'species', 'vital_rate', 'lambda')) %>% 
  group_by(species, vital_rate, model) %>% 
  arrange( vital_rate, species, model , lambda) %>% 
  filter( n() > 1) %>% 
  mutate( mean_lppd = mean(lppd, na.rm = TRUE), 
          sd_lppd = sd(lppd, na.rm = TRUE), 
          outlier = ifelse( is.na(lppd) | abs(lppd - mean_lppd) > 400 , TRUE, FALSE ), 
          max_lppd = max(lppd[!outlier])) %>% 
  group_by( species, vital_rate, model, outlier) %>% 
  mutate( lfit = ifelse( outlier, NA, loess_predict(lppd , lambda )), 
          lfit.se = ifelse( outlier, NA, loess_predict(lppd, lambda, TRUE)), 
          lppd_diff = lppd - max_lppd )

gg_lppd_base <- 
  ggplot( lppds, aes( x = lambda, y = lppd )) + 
  geom_text( aes( x = lambda, y = lppd, label = lambda), size = 3) + 
  geom_point( data = . %>% filter( lppd == max(lppd)), aes( x = lambda, y = lppd ), size = 4.7, color = 'red', shape = 1) 

gg_no_outlier <- 
  lppds %>% 
  filter (!outlier ) %>% 
  group_by(species, vital_rate, model) %>% 
  filter( n() > 1 ) %>% 
  do(gg = gg_lppd_base %+% . + 
       geom_smooth(aes( group = 1), color = 'blue',  se = FALSE) + 
       ggtitle(paste('out of sample lppd by regularization plot for', unique(.$species), unique(.$vital_rate), 'model', unique(.$model)) ))

gg <- lppds %>% 
  group_by(species, vital_rate, model) %>% 
  filter( n() > 1 ) %>% 
  do(gg = gg_lppd_base %+% . + 
       ggtitle(paste('out of sample lppd by regularization plot for', unique(.$species), unique(.$vital_rate), 'model', unique(.$model)) ))

# output ----------------------------------------------------------------------------------------

pdf('figures/plot_lppd_by_lambda_no_outliers.pdf', height = 8 , width = 8 )

print( gg_no_outlier$gg ) 

dev.off()


pdf('figures/plot_lppd_by_lambda.pdf', height = 8 , width = 8 )

print( gg$gg ) 

dev.off()

# save highest lppd models 

best_fits <- 
  lppds %>% 
  group_by( vital_rate, species, model ) %>% 
  filter( lfit == max(lfit[!is.na(lfit)] )) 

best_fits <- left_join( best_fits, model_table)

best_fits <- best_fits %>% arrange(vital_rate, species, model )

write.csv(best_fits , 'output/best_lppd_scores.csv')
