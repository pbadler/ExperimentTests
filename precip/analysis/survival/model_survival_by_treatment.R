rm(list = ls())
library(lme4)

model_survival_by_treatment <- function( allD , species_name, output_dir) { 
  
  allD <- subset(allD, Treatment %in% c('Control', 'Drought', 'Irrigation') )
  allD <- subset(allD , year > 2006)
  
  allD$year <- as.factor(allD$year)

  # use lmer
  # w/o treatment effect
  
  m0 <- glmer(survives~logarea + W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts +
               (1|Group)+(logarea|year),data=allD, family = 'binomial') 
  
  # w/ treatment effect
  m1 <- glmer(survives~logarea+Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts +
               (1|Group)+(logarea|year),data=allD, family = 'binomial') 
  # 
  print( paste('RESULTS FOR:', species_name ))
  print( anova(m1, m0) )  # Test treatment effect 
  
  lmer_results <- list(m0, m1)
  
  # save output -------------------------------------------------------------------
  
  saveRDS(lmer_results, file = file.path(output_dir, paste0(species_name, 'survival_treatment_effects.lmer.RDS')))
  
  # -------------------------------------------------------------------------------
} 


#########################################
#  2. Fit survival models 
#########################################

survival_data <- dir('data/temp_data/', pattern = 'survival.RDS', full.names = TRUE)

spp_names <- regmatches(survival_data, m = gregexpr( pattern = '([A-Z]{4})', survival_data )) 

allD <- lapply(survival_data, readRDS)

output_dir <- 'output/'


for( i in 1:length(survival_data)){ 
  
  model_survival_by_treatment(allD = allD[[i]], spp_names[[i]], output_dir = output_dir)
  
} 
