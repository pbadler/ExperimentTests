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

load('data/temp_data/master_list.Rdata')

lppd_scores <- 
  all_lppds %>% 
  filter( type == 'out_of_sample' ) %>%
  mutate( clean_fn = str_replace(fn, pattern = '_WAIC\\.csv$', replacement = '')) %>%
  separate(clean_fn , into = c('species', 'vital_rate', 'model', 'lambda', 'chains', 'year_oos'), sep = '_')  

model_table <- model_table %>% filter( vital_rate == 'growth')

lppd_scores <- merge(model_table, lppd_scores , by = c('model' , 'species', 'vital_rate', 'lambda', 'year_oos'), all.x = TRUE) 

leave_one_out <- 
  lppd_scores %>% 
  filter( year_oos != 'c(2007:2015)') %>% 
  arrange( species, vital_rate, lambda, year_oos ) %>% 
  select( lppd, model, species, vital_rate, lambda, year_oos) %>% 
  group_by( species, vital_rate, model, lambda ) %>% 
  summarise(lppd = mean(lppd)) 

true_hold_out <- 
  lppd_scores %>% 
  filter( year_oos == 'c(2007:2015)') %>% 
  group_by( species, vital_rate, model, lambda ) %>% 
  summarise(lppd = sum(lppd) , n = n() )

leave_one_out <- 
  leave_one_out %>% 
  group_by(species, vital_rate, model) %>% 
  arrange( vital_rate, species, model , lambda) 


lo <- split(leave_one_out, paste0(leave_one_out$species, leave_one_out$vital_rate))

for(i in 1:length(lo)){
  plot(data = lo[[i]], lppd ~ lambda) 
}

plot(data = leave_one_out %>% filter( species == 'PSSP', vital_rate == 'growth'), lppd ~ lambda)
plot(data = lppd_scores %>% filter(species == 'PSSP', vital_rate == 'growth' , year_oos != 'c(2007:2015)'), lppd ~ lambda)

lppd_scores %>% filter( species == 'PSSP', vital_rate == 'growth', year_oos != 'c(2007:2015)' , lppd < -300) %>% group_by(lambda) %>% summarise( mlppd = mean(lppd), n = n())

test <- lppd_scores %>% filter (species == 'POSE', vital_rate == 'growth', year_oos != 'c(2007:2015)')

test <- readRDS('data/temp_data/PSSP_growth.RDS')
table(test$year, test$Grazing)

# ----------------------------------------------------------------------------------------------------------------# 

gg_lppd_base <- 
  ggplot( leave_one_out, aes( x = lambda, y = lppd )) + 
  geom_text( aes( x = lambda, y = lppd, label = lambda))  

gg_lppd <- 
  gg_lppd_base +  
  geom_line(aes(x = lambda, y = lfit - max_lppd, group = 1))


# gg_no_outlier <- 
#   lppds %>% 
#   filter (!outlier ) %>% 
#   group_by(species, vital_rate, model) %>% 
#   filter( n() > 1 ) %>% 
#   do(gg = gg_lppd_base %+% . + 
#        geom_smooth(aes( group = 1), color = 'blue',  se = FALSE) + 
#        ggtitle(paste('out of sample lppd by regularization plot for', unique(.$species), unique(.$vital_rate), 'model', unique(.$model)) ))

gg <- 
  leave_one_out %>% 
  group_by(species, vital_rate, model) %>% 
  filter( n() > 1 ) %>% 
  do(gg = gg_lppd_base %+% . + 
       ggtitle(paste('leave one out, lppd by regularization plot for', unique(.$species), unique(.$vital_rate), 'model', unique(.$model)) ))

# output ----------------------------------------------------------------------------------------

dfs <- split( leave_one_out, paste0(leave_one_out$species, leave_one_out$vital_rate) )
length(dfs)
dfs
dfs[[1]]
plot( dfs[[1]]$lambda , dfs[[1]]$lppd ) 


pdf('figures/plot_lppd_by_lambda_leave_one_out.pdf', height = 8 , width = 8 )

print( gg$gg ) 

dev.off()

# ----------------------------------------------------------------------------------------------- 

gg_out <- 
  true_hold_out %>% 
  group_by(species, vital_rate) %>% 
  filter( n() > 1 ) %>% 
  do(gg = gg_lppd %+% . + 
       ggtitle(paste('out of sample lppd by regularization plot for', unique(.$species), unique(.$vital_rate), 'model', unique(.$model)) ))

# save highest lppd models 

best_fits <- 
  leave_one_out %>% 
  group_by( vital_rate, species, model ) %>% 
  filter( lfit == max(lfit[!is.na(lfit)] )) 

no_reg_fits <- 
  all_lppds %>% 
  filter( model == 1 ) 

best_fits <- bind_rows ( best_fits %>% select_(.dots = names( all_lppds)), no_reg_fits )
best_fits <- left_join( best_fits, model_table)

best_fits <- best_fits %>% arrange(vital_rate, species, model )

write.csv(best_fits , 'output/best_lppd_scores_oos.csv')
