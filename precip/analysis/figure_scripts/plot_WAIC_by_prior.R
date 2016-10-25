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
  function( waic , lambda, error = FALSE ) { 
    test_fit <- loess(formula = waic ~ lambda)
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

# ----------------------------------------------------------------------------------------------------------------# 

all_waics <- 
  all_waics %>% 
  mutate( clean_fn = str_replace(fn, pattern = '_WAIC\\.csv$', replacement = '')) %>%
  separate(clean_fn , into = c('species', 'vital_rate', 'model', 'lambda', 'chains'), sep = '_') %>% 
  mutate( model = as.numeric(model), lambda = as.numeric( lambda))

waics <- 
  merge(all_waics, model_table, by = c('model', 'species', 'vital_rate', 'lambda')) %>% 
  group_by(species, vital_rate, model) %>% 
  arrange( vital_rate, species, model , lambda) %>% 
  filter( n() > 1 ) %>% 
  mutate( mean_WAIC = mean(waic), 
          sd_WAIC = sd(waic), 
          outlier = ifelse(abs(waic - mean_WAIC) > 400 , TRUE, FALSE ), 
          min_waic = min(waic[!outlier])) %>% 
  group_by( species, vital_rate, model, outlier) %>% 
  mutate( lfit = ifelse( outlier, NA, loess_predict(waic , lambda )), 
          lfit.se = ifelse( outlier, NA, loess_predict(waic, lambda, TRUE)), 
          waic_diff = waic - min_waic )

gg_waic_base <- 
  ggplot( waics, aes( x = lambda, y = waic_diff )) + 
  geom_text( aes( x = lambda, y = waic_diff, label = lambda))  

gg_waic <- 
  gg_waic_base +  
  geom_line(aes(x = lambda, y = lfit - min_waic, group = 1))

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

write.csv(outliers, 'output/outlier_runs.csv', row.names = FALSE)

# save lowest WAIC models 

waics$lfit

best_fits <- 
  waics %>% 
  group_by( vital_rate, species, model ) %>% 
  filter( lfit == min(lfit[!is.na(lfit)] )) 

no_reg_fits <- 
  all_waics %>% 
  filter( model == 1 ) 

best_fits <- bind_rows ( best_fits %>% select_(.dots = names( all_waics)), no_reg_fits )
best_fits <- left_join( best_fits, model_table)

best_fits <- best_fits %>% arrange(vital_rate, species, model )

write.csv(best_fits , 'output/best_WAIC_scores.csv')