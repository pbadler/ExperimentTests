###########################################################################################
#
# Read files list  
# gather WAIC scores in dataframe 
# save 
#
###########################################################################################

rm(list = ls())


# input ------------------------------------------------------------------------------------
files <- dir('../output/WAIC_scores', pattern = '*WAIC.csv', full.names = T)

# collect ----------------------------------------------------------------------------------

WAIC_scores <- do.call( rbind , lapply( files, read.csv)  ) 

# output ------------------------------------------------------------------------------------

write.csv(WAIC_scores , '../output/WAIC_scores_compiled.csv', row.names = FALSE)
