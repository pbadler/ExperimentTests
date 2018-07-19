######################################################################################
#
# Make STAN datalist  
#
#####################################################################################

# These retrieve and aggregate all the climate data 

#source('R/climate/ExtractData_3Runs.R') # don't run this because using the deprecated version of the soilwat package required is a pain in the ass 
source('R/climate/aggregate_spot_VWC.R')
source('R/climate/soilMoistureTreatmentEffects.R')
source('R/climate/aggregate_VWC_data.R')

# ----------------------------------------------------------------- 
source('R/get_all_demographic_data.R')
source('R/climate/make_climate_variables.R')
source('R/climate/prepare_climate_covariates.R')

source('R/prep_vital_rate_df.R')

