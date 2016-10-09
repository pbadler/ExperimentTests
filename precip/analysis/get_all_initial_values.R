#

# get all initial values 
rm(list = ls() )

scripts <- dir('analysis', 'get_initial_values.R', recursive = T, full.names = T)


lapply( scripts, source)


