############################################################# 
#
# NOTE: The whole thing takes several hours to run. (most of day?)
# Involves three rounds of model fitting. 
#
#############################################################

rm(list = ls() ) 

setwd("C:/Repos/ExperimentTests/precip/")
my_path <- "c:/repos/" # this is needed to get to the driversdata directory

library(rstan)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(lme4)
library(zoo)
library(texreg)
library(xtable)
library(gridExtra)
library(MASS)

# run data preparation files first --------------------------- # 

source('analysis/figure_scripts/save_plot_theme.R')

source('R/prepare_datalists_for_STAN.R')

rm(list=setdiff(ls(), "my_path")) # clean up, but leave my_path

# analysis pipeline ------------------------------------------ # 

# 1. Soil moisture analysis 
source('analysis/figure_scripts/plot_spring_soil_moisture_spot_measures.R')
rm(list=setdiff(ls(), "my_path")) # clean up, but leave my_path

# 2. Cover trends

source('analysis/treatment_trends_precip.R')
rm(list=setdiff(ls(), "my_path")) # clean up, but leave my_path

# 3. fit treatment models 

source('analysis/fit_stan_treatment_models.R') # takes a few hours to run 

source('analysis/figure_scripts/plot_treatment_stan_fits.R') # takes a few minutes to run 

source('analysis/summarize_treatment_model_fits.R') # take several minutes

source('analysis/figure_scripts/plot_treatment_parameter_estimates.R')

# 4. run year effects model 
  
source('analysis/fit_stan_year_effects_models.R')

source('analysis/summarize_year_model_fits.R')
source('analysis/figure_scripts/plot_year_model_fits.R')

# 5. find climate correlations 

source('analysis/find_climate_year_effects_correlations.R')

source('analysis/select_climate_covariates.R')


# 6. modify datalists with selected climate covariates

source('analysis/modify_climate_datalists.R')

# 7. fit climate models 

source('analysis/fit_stan_climate_models.R')

source('analysis/check_for_divergent_transitions.R')

source('analysis/figure_scripts/plot_climate_fits.R')

source('analysis/figure_scripts/plot_compare_year_model_climate_model_year_effects.R')

source('analysis/summarize_climate_model_fits.R')

source('analysis/figure_scripts/plot_climate_parameter_estimates.R')

# 8. get prediction scores 

source('analysis/calculate_lppd.R')
source('analysis/calculate_total_lppd.R')
source('analysis/get_MSE_scores.R')

source('analysis/make_prediction_score_tables.R')

# 9. plot observed and predicted climate effects on vital rates 

source('analysis/figure_scripts/plot_obs_vs_pred_treatment_parameters.R')

source('analysis/figure_scripts/plot_significant_effect_comparison.R')

# 8. plot observed cover and predicted cover 

source('analysis/figure_scripts/plot_one_step_ahead_cover.R')

source('analysis/figure_scripts/plot_predicted_vs_observed_cover_changes.R')

source('analysis/figure_scripts/plot_cover.R')


# 9. 

