####################################################################################
#
# run scripts to fetch all demographic data 
# 
####################################################################################

get_growth_files <- dir('R', pattern = 'get.*growth\\.R', recursive = TRUE, full.names = TRUE)

get_survival_files <- dir('R', pattern = 'get.*survival\\.R', recursive = TRUE, full.names = TRUE)

get_recruitment_files <- dir('R', pattern = 'get.*recruitment\\.R', recursive = TRUE, full.names = TRUE)

lapply( c(get_growth_files, get_survival_files, get_recruitment_files), source) 

#TODO write functions to export to HPC server 

