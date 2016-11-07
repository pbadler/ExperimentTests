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


all_lppds <- read.csv('output/WAIC_scores_oos.csv')

model_table <- read.csv('data/temp_data/model_table_oos.csv')

lppd_scores <- 
  all_lppds %>% 
  filter( type == 'out_of_sample') %>%
  mutate( clean_fn = str_replace(fn, pattern = '_WAIC\\.csv$', replacement = '')) %>%
  separate(clean_fn , into = c('species', 'vital_rate', 'model', 'lambda', 'chains', 'year_oos'), sep = '_')  

lppd_scores <- merge(model_table, lppd_scores , by = c('model' , 'species', 'vital_rate', 'lambda', 'year_oos'), all.x = TRUE) 

lppd_scores <- 
  lppd_scores %>% 
  filter( year_oos != 'c(2007:2015)') %>% 
  group_by( species, vital_rate, model, lambda ) %>% 
  summarise(lppd = sum(lppd) , n = n() )

load('data/temp_data/master_list.Rdata')
# ----------------------------------------------------------------------------------------------------------------# 

# lppds <- 
#   merge(all_lppds, model_table, by = c('model', 'species', 'vital_rate', 'lambda')) %>% 
#   group_by(species, vital_rate, model) %>% 
#   arrange( vital_rate, species, model , lambda) %>% 
#   filter( n() > 1) %>% 
#   mutate( mean_lppd = mean(lppd, na.rm = TRUE), 
#           sd_lppd = sd(lppd, na.rm = TRUE), 
#           outlier = ifelse( is.na(lppd) | abs(lppd - mean_lppd) > 400 , TRUE, FALSE ), 
#           max_lppd = max(lppd[!outlier])) %>% 
#   group_by( species, vital_rate, model, outlier) %>% 
#   mutate( lfit = ifelse( outlier, NA, loess_predict(lppd , lambda )), 
#           lfit.se = ifelse( outlier, NA, loess_predict(lppd, lambda, TRUE)), 
#           lppd_diff = lppd - max_lppd )

gg_lppd_base <- 
  ggplot( lppd_scores, aes( x = lambda, y = lppd )) + 
  geom_text( aes( x = lambda, y = lppd, label = lambda))  

gg_lppd <- 
  gg_lppd_base +  
  geom_smooth(aes(x = lambda, y = lppd, group = 1), se = FALSE)

# gg_no_outlier <- 
#   lppd_scores %>% 
#   filter (!outlier ) %>% 
#   group_by(species, vital_rate, model) %>% 
#   filter( n() > 1 ) %>% 
#   do(gg = gg_lppd_base %+% . + 
#        geom_smooth(aes( group = 1), color = 'blue',  se = FALSE) + 
#        ggtitle(paste('out of sample lppd by regularization plot for', unique(.$species), unique(.$vital_rate), 'model', unique(.$model)) ))

gg <- 
  lppd_scores %>% 
  group_by(species, vital_rate) %>% 
  filter( n() > 1 ) %>% 
  do(gg = gg_lppd %+% . + 
       ggtitle(paste('out of sample lppd by regularization plot for', unique(.$species), unique(.$vital_rate), 'model', unique(.$model)) ))

#gg_lppd %+% (lppd_scores %>% filter( species == 'ARTR', vital_rate == 'growth'))

gg$gg[[12]]

#gg_lppd %+% subset( lppd_scores, species == 'ARTR'  & vital_rate == 'growth')

#dev.off()


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

no_reg_fits <- 
  all_lppds %>% 
  filter( model == 1 ) 

best_fits <- bind_rows ( best_fits %>% select_(.dots = names( all_lppds)), no_reg_fits )
best_fits <- left_join( best_fits, model_table)

best_fits <- best_fits %>% arrange(vital_rate, species, model )

write.csv(best_fits , 'output/best_lppd_scores.csv')
