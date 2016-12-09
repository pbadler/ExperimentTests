# 
#
# construct model table 
#
#
source('analysis/make_master_model_parameter_list.R')

rm(list = ls())

library(dplyr)
library(tidyr)

load('data/temp_data/master_list.Rdata')

print( 'model parameter list:')
print(master_list)

# survival model sets with parameters  ----------------------------------------- 

models <- expand.grid( 
  species = master_list$species,
  vital_rate = master_list$vital_rates$vital_rate, 
  model = master_list$smodels$model, 
  lambda = c(1:master_list$nlambda))

all_pars <- rbind( master_list$smodels, master_list$gmodels, master_list$rmodels)

models <- merge(models, all_pars, by = c('vital_rate', 'model'))

models <- merge( models, master_list$sd_vec)
models <- merge( models, master_list$vital_rates)

models <- models %>% arrange( vital_rate, species, model, lambda )

full_table <- 
  models %>% 
  arrange( vital_rate, species, model, lambda ) %>% 
  mutate( nlambda = master_list$nlambda, 
          index = row_number())

print (paste0( 'table has ' , nrow( full_table ), ' rows,   does this match expected??? '))

short_table <- 
  full_table %>% 
  filter( lambda == 1 )

write.csv( full_table, 'data/temp_data/model_table.csv', row.names = FALSE)

write.csv( short_table, 'data/temp_data/short_model_table.csv', row.names = FALSE )


# spp_tables <- split(full_table, interaction(full_table$species, full_table$vital_rate, full_table$lambda ))
# 
# for( i in 1:length(spp_tables)){ 
#   
#   species <- spp_tables[[i]]$species
#   vr <- spp_tables[[i]]$vital_rate
#   
#   write.csv(paste0 ( 'data/temp_data/regularization_tables_', species, '_', vital_rate, '_table.csv' ) )
# }


