# 
#
# construct model table 
#
#
library(dplyr)
library(tidyr)

rm(list = ls())

do_spp <- c('ARTR','HECO', 'POSE', 'PSSP')
do_vr_number <- c(1, 2, 3)
do_vr <-  c('growth', 'survival', 'recruitment')
do_model <- c(1:5)
do_parameters <- c()

do_prior <- c(1:30)

models <- expand.grid(species = do_spp, vital_rate_number = do_vr_number, model = do_model, prior = do_prior)

nrow( models ) 

vital_rates <- data.frame(vital_rate = do_vr, vital_rate_number = do_vr_number)

models <- merge(models, vital_rates)

models <- 
  models %>% 
  arrange( vital_rate_number, species, model, prior ) %>% 
  filter( !(model == 1 & prior > 1) ) %>%
  filter( !(model == 3 & prior > 1) )


write.csv( models, 'data/temp_data/model_table.csv', row.names = FALSE)


