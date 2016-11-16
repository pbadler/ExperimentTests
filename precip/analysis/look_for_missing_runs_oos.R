#
# Find model combinations that did not run 
#
#

library(dplyr)
library(tidyr)
library(stringr)

waics <- read.csv('output/WAIC_scores_oos.csv')

model_table <- read.csv('data/temp_data/model_table_oos.csv')

waics <- 
  waics %>% 
  filter( type == 'out_of_sample') %>%
  mutate( clean_fn = str_replace(fn, pattern = '_WAIC\\.csv$', replacement = '')) %>%
  separate(clean_fn , into = c('species', 'vital_rate', 'model', 'lambda', 'chains', 'year_oos'), sep = '_')  

models_missing <- 
    merge(model_table, waics , by = c('model' , 'species', 'vital_rate', 'lambda', 'year_oos'), all.x = TRUE) %>% filter( is.na(lppd))

models_missing$d <- 1:nrow(models_missing)

models_missing <- split( models_missing, ceiling( seq_along(models_missing$d)/900)) # split into  chunks of 900 

for( i in 1:length(models_missing)) { 

  write.csv(models_missing[[i]] %>% 
        arrange(index), file = paste0( 'output/missing_runs_oos', i, '.csv'), row.names = FALSE) 
}


