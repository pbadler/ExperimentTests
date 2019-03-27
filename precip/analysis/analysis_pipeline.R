############################################################# 
#
# NOTE: Running all the K-fold cross validation takes days 
#
#############################################################

rm(list = ls() ) 

setwd("~/Dropbox/projects/ExperimentTests/precip")
my_path <- "~" # this is needed to get to the driversdata directory

library(rstan)
library(tidyverse)
library(stringr)
library(zoo)
library(texreg)
library(xtable)
library(gridExtra)
library(MASS)
library(ggmcmc)

# run data preparation files first --------------------------- # 

source('analysis/figure_scripts/save_plot_theme.R')

source('R/prepare_data_for_STAN.R')

rm(list=setdiff(ls(), "my_path")) # clean up, but leave my_path

# analysis pipeline ------------------------------------------ # 

# 1. Soil moisture analysis 
source('analysis/figure_scripts/plot_spring_soil_moisture_spot_measures.R')
rm(list=setdiff(ls(), "my_path")) # clean up, but leave my_path

# 2. Cover trends

source('analysis/treatment_trends_precip.R')
rm(list=setdiff(ls(), "my_path")) # clean up, but leave my_path

# ----- Model Fitting/Selection ---------------------------------- # 

# 3. Fit candidate models and evaluate using K-Fold cross validation

source('analysis/fit_growth_models1.R') # takes ~ a while 

source('analysis/fit_survival_models1.R') # takes ~ a while 

source('analysis/fit_recruitment_models1.R') # takes ~ a while 

# 4. Rank models based on LPPD 

source('analysis/rank_models.R')
source('analysis/annotate_model_ranks.R')

# 5. Re-run top models for each species and vital rate 


# ------ Simulation ---------------------------------------------  #

# 6. Generate IBM predictions based on top demographic models 


# ----- Validation ----------------------------------------------- # 

# 7. Evaluate vital rate predictions on the held-out validation data (2012 to 2016)

# 8. Evalute cover predictions from IBMs 


# ----- Generate Figures ----------------------------------------- # 



# ----- Generate Tables  ----------------------------------------- # 


# ----- Knit manuscript  ------------------------------------------#

