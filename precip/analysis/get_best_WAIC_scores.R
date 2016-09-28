# 
# Find best prior based on WAIC 
# Save output as table that can be run by the run stan models simple script  
# 

rm(list = ls())

df <- read.csv('output/WAIC_scores.csv')

model_table <- read.csv('data/temp_data/short_model_table.csv')

df$clean_fn <- gsub(as.character(df$fn), pattern = '_WAIC\\.RDS$', replacement = '')

id <- data.frame( do.call( rbind, strsplit(as.character(df$clean_fn), split = '_') ) )  
names(id ) <- c('species', 'vital_rate', 'model', 'prior', 'chains')

df <- data.frame( df, id )

df_list <- split( df , factor(paste(df$species, df$vital_rate, df$model)))

best_waic <- do.call( rbind, lapply( df_list, function(x) x[which.min(x$waic), ]))

best_waic$vital_rate <- as.character( best_waic$vital_rate ) 
model_table$vital_rate <- as.character(model_table$vital_rate)

model_table <- model_table[ , - grep (names(model_table), pattern = 'prior') ] # drop old prior column 

best_waic <- merge( best_waic, model_table, all.x = TRUE, by = c('species', 'model', 'vital_rate'))

best_waic <- best_waic[ order(best_waic$vital_rate, best_waic$species, best_waic$model), ]

write.csv(best_waic, 'output/best_WAIC_scores.csv', row.names = FALSE)


