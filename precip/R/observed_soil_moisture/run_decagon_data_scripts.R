###############################################################################
# 
# Run all scripts related to decagon data processing 
# 
###############################################################################

home <- '~'
# setwd(file.path(home, 'driversdata', 'data', 'idaho_modern', 'soil_moisture'))
setwd(file.path( home, 'Documents','ExperimentTests', 'precip', 'R'))

source('import_and_format_decagon_data.R')

source('R/observed_soil_moisture/correct_decagon_dates.R')

source('R/observed_soil_moisture/correct_decagon_readings.R')

source('R/observed_soil_moisture/merge_decagon_with_climate_station_data.R')

source('R/observed_soil_moisture/export_daily_soil_moisture_for_SOILWAT.R')
