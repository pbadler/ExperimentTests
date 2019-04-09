######################################################################################
#
# Make STAN datalist  
#
#####################################################################################

# These retrieve and aggregate all the climate data 

# source('R/climate/ExtractData_3Runs.R') # depends on soilwat data sent by Caitlin Andrews

# don't run this unless you can install Rsoilwat31, an old version of
# the soilwat package. 
# Instructions for installing the version are given in the script.  


source('R/climate/aggregate_spot_VWC.R')
source('R/climate/soilMoistureTreatmentEffects.R')
source('R/climate/aggregate_VWC_data.R')

# ----------------------------------------------------------------- 
source('R/get_all_demographic_data.R') # depends on access to driversdata 
source('R/calculate_cover_per_plot.R') # depends on access to driversdata 

source('R/climate/make_climate_variables.R') 
source('R/climate/prepare_climate_covariates.R')

source('R/prep_vital_rate_df.R')
