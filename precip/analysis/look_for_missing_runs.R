#
# Find model combinations that did not run 
#
#

library(dplyr)
library(tidyr)
library(stringr)

waics <- read.csv('output/WAIC_scores.csv')

model_table <- read.csv('data/temp_data/model_table.csv')

waics <- 
  waics %>% 
  mutate( clean_fn = str_replace(fn, pattern = '_WAIC\\.csv$', replacement = '')) %>%
  separate(clean_fn , into = c('species', 'vital_rate', 'model', 'prior', 'chains'), sep = '_')  

models_missing <- merge(model_table, waics , by = c('model' , 'species', 'vital_rate', 'prior'), all.x = TRUE)


print( models_missing %>% 
  filter( vital_rate != 'recruitment') %>% 
  filter( is.na(waic) ) %>% 
  arrange(index) ) 

