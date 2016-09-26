# 
# Find best prior based on WAIC 
# 
# 

rm(list = ls())

df <- readRDS('../output/WAIC_scores.RDS')

df$clean_fn <- gsub(as.character(df$fn), pattern = '^\\./', replacement = '')
df$clean_fn <- gsub(as.character(df$clean_fn), pattern = '\\.RDS$', replacement = '')

id <- data.frame( do.call( rbind, strsplit(as.character(df$clean_fn), split = '_') ) )  
names(id ) <- c('species', 'vital_rate', 'model', 'prior', 'chains')

df <- data.frame( df, id )

df_list <- split( df , factor(paste(df$species, df$vital_rate, df$model)))

best_waic <- do.call( rbind, lapply( df_list, function(x) x[which.min(x$waic), ]))

saveRDS(best_waic, '../output/best_WAIC_scores.RDS')


