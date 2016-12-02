# run data preparation files first --------------------------- # 

# analysis pipeline ------------------------------------------ # 

# 1. run year effects model 
  
# source('analysis/fit_year_effects.R')
# source('analysis/summarize_year_model_fits.R')

# 2. find climate correlations 

# source('analysis/find_climate_year_effects_correlations.R')

# source('analysis/select_climate_covariates.R')

# 3. modify datalists with selected climate covariates

# source('analysis/modify_climate_datalists.R')

# 4. fit climate models 

# source('analysis/fit_stan_climate_models.R')

# source('analysis/summarize_climate_model_fits.R')

# source('analysis/figure_scripts/plot_predictive_fits.R')

# 5. fit treatment models 

# source('analysis/fit_treatment_effects.R')

# source('analysis/summarize_treatment_model_fits.R')

# source('analysis/figure_scripts/plot_treatment_stan_fits.R')

# source('analysis/figure_scripts/plot_treatmeant_parameter_estimates.R')
 
# 6.  check for divergence 

# source('analysis/check_for_divergent_transitions.R')

# 7. plot observed cover and predicted cover 

# 8. plot observed and predicted climate effects on vital rates 

# 9. 