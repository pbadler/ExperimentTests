# 
#
# construct model table 
#
#
rm(list = ls())

library(dplyr)
library(tidyr)

load('data/temp_data/master_list.Rdata')

print( 'model parameter list:')
print(master_list)

models <- expand.grid(
  species = master_list$species, 
  vital_rate = master_list$vital_rates$vital_rate, 
  model = master_list$models$model, 
  lambda = c(1:master_list$nlambda))

models <- merge(models, master_list$models)

models <- merge( models, master_list$sd_vec)

models <- merge( models, master_list$vital_rates)

models <- models %>% arrange( vital_rate, species, model, lambda )

rmodels <- expand.grid( 
  species = master_list$species, 
  vital_rate = 'recruitment',
  model = master_list$models$model, 
  lambda = c(1:master_list$nlambda))

rmodels <- merge( rmodels, master_list$rmodels)

rmodels <- merge( rmodels, master_list$sd_vec)

rmodels <- merge( rmodels, master_list$vital_rates)

models <- models[ !models$vital_rate == 'recruitment', ]

models <- rbind( models, rmodels)

models <- models %>% 
  arrange( vital_rate, species, model, lambda)

full_table <- 
  models %>% 
  arrange( vital_rate, species, model, lambda ) %>% 
  filter( !(model == 1 & lambda != ceiling(0.5*master_list$nlambda) )) %>% 
  mutate( nlambda = master_list$nlambda, 
          index = row_number())

print (paste0( 'table has ' , nrow( full_table ), ' rows,   does this match expected??? '))

short_table <- 
  full_table %>% 
  filter( lambda == ceiling(0.5*master_list$nlambda) )

write.csv( full_table, 'data/temp_data/model_table.csv', row.names = FALSE)

write.csv( short_table, 'data/temp_data/short_model_table.csv', row.names = FALSE )



