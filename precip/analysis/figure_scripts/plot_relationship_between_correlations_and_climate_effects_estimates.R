# 
library(stringr)
rm(list = ls())

rcor <- read.csv('output/year_effects_correlations_recruitment.csv')
gscor <- read.csv('output/year_effects_correlations_growth_survival.csv')

clim_model_fs <- dir('output', 'climate_model_parameter.*.csv', full.names = T)

covariates <- read.csv('output/selected_climate_covariates.csv')

all_clim_model_pars <- lapply( clim_model_fs, read.csv)

all_clim_model_pars <- do.call( rbind, all_clim_model_pars)

all_clim_model_pars <- all_clim_model_pars[ str_detect ( all_clim_model_pars$X, 'b2' ), ]

all_clim_model_pars <- 
  all_clim_model_pars %>% 
  arrange( species, vital_rate)

covariates <- 
  covariates %>% 
  arrange(species, vital_rate)

covariate_names <- str_split(covariates$covars, ',')

all_clim_model_pars$parameter <- unlist(covariate_names)

# get correlations with year effects 
gscor <- 
  gscor %>% 
  gather( correlation, val , Intercept, size ) %>%
  select( -Int_pval, -size_pval)

rcor <- 
  rcor %>% 
  gather(correlation, val, Intercept ) %>% 
  select( -Int_pval) 

head(gscor)
head(rcor)
vrcor <- rbind(gscor, rcor)

vrcor$parameter <- as.character(vrcor$X )

vrcor <- 
  vrcor %>% 
  mutate( parameter = ifelse(correlation == 'size', paste0(parameter, 'xlogarea.t0'), parameter ) )


all_vr_pars <- merge(vrcor, all_clim_model_pars, by = c('species', 'vital_rate', 'parameter'))

plot( all_vr_pars$val, all_vr_pars$mean ) # plot relationship between correlation with year effects and climate parameter estimate
points(data = subset(all_vr_pars, correlation == 'size'), mean ~ val , col = 'red')
abline(0,1, lty = 2)

plot( all_vr_pars$val, all_vr_pars$mean ) # plot relationship between correlation with year effects and climate parameter estimate
points(data = subset(all_vr_pars, vital_rate == 'growth'), mean ~ val , col = 'red')
abline(0,1, lty = 2)
