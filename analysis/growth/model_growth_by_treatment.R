rm(list = ls())
library(lme4)

model_growth_by_treatment <- function( allD , species_name, output_dir) { 
  
  allD <- subset(allD, Treatment %in% c('Control', 'Drought', 'Irrigation') )
  allD <- subset(allD , year > 2006)
  
  allD$year <- as.factor(allD$year)

  # use lmer
  # w/o treatment effect
  
  m0 <- lmer(logarea.t1~logarea.t0+W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts +
               (1|Group)+(logarea.t0|year),data=allD) 
  
  # w/ treatment effect
  m1 <- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts +
               (1|Group)+(logarea.t0|year),data=allD) 
  # 
  print( paste('RESULTS FOR:', species_name ))
  print( anova(m1, m0) )  # Test treatment effect 
  
  lmer_results <- list(m0, m1)
  
  # save output -------------------------------------------------------------------
  
  saveRDS(lmer_results, file = file.path(output_dir, paste0(species_name, 'growth_treatment_effects.lmer.RDS')))
  
  # -------------------------------------------------------------------------------
} 


#########################################
#  2. Fit growth models 
#########################################

growth_data <- dir('data/temp_data/', pattern = 'growth.RDS', full.names = TRUE)

spp_names <- regmatches(growth_data, m = gregexpr( pattern = '([A-Z]{4})', growth_data )) 

allD <- lapply( growth_data, readRDS)

output_dir <- 'output/'

for( i in 1:length(growth_data)){ 
  
  model_growth_by_treatment(allD = allD[[i]], spp_names[[i]], output_dir = output_dir)
  
} 
