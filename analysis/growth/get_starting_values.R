##################################################################################################################
#
# Fit models to historical data only 
#
# Export starting values
#
##################################################################################################################

rm(list = ls())

library(lme4)

dl <- readRDS('data/temp_data/all_growth_combined.RDS')

test <- dl[[1]]

# climate effects 

test$P.1.2 # year two, spring precip  
test$P.4.2 # year two, total annual precip  
test$P.4.1 # year one, total annual precip 

test$T.1.2 # year two, spring temp  
test$T.1.1 # year one, spring temp 
test$T. 

#---model descriptions --------------------------------------------------------------------------------------------
#  
#   m0: null model, w/o climate, w/o competition 
#   m1: climate model, w/ climate, w/o competition 
#   m2: single species model, w/o climate, w/ intra-specific competition 
#   m3: single species climate model, w/ climate, w/ intra-specific competition
#   m4: community model, w/climate, w/ intra and interspecific competition 
#
# ----------------------------------------------------------------------------------------------------------------- 

m0 <- lmer( logarea.t1 ~ logarea.t0 + (1|Group) + (logarea.t0|year), data = test ) 

m1 <- lmer( logarea.t1 ~ logarea.t0 + 
              P.1.1 + # current 
              
              Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts +
             (1|Group)+(logarea.t0|year),data=allD) 

m1 <- lmer(logarea.t1 ~ logarea.t0 + W.ARTR + W.HECO + W.POSE + W.PSSP + W.allcov + W.allpts +  # competition effects
             (1|Group) + 
             (logarea.t0|year), data=allD) 

m1 <- lmer(logarea.t1~logarea.t0+Treatment+W.ARTR + W.HECO + W.POSE + W.PSSP+  W.allcov + W.allpts +
             (1|Group)+(logarea.t0|year),data=allD) 