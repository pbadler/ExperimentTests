#
# test for period effects 
# 
rm(list = ls())
library(lme4)
library(tidyr)
library(dplyr)
clim <- readRDS('data/temp_data/all_clim_covs.RDS')
clim <- subset(clim, year > 1926)

clim[ , grep( '^[PT]\\.|^VWC', names(clim))] <- scale( clim[ , grep( '^[PT]\\.|^VWC', names(clim))])

# ag_VWC <- 
#   clim %>% 
#   gather( var, val, starts_with( 'VWC')) %>% 
#   separate(var, c('var', 'layer'), '_') %>% 
#   group_by(Period, Treatment, year ,quarter, var) %>%  
#   summarise( val = mean(val)) %>% 
#   spread( var, val )

# clim <- merge( clim %>% select( - starts_with('VWC')), ag_VWC )

ARTR <- readRDS('data/temp_data/ARTR_survival.RDS')
HECO <- readRDS('data/temp_data/HECO_survival.RDS')
POSE <- readRDS('data/temp_data/POSE_survival.RDS')
PSSP <- readRDS('data/temp_data/PSSP_survival.RDS')

# 
pairs( clim[, grep('sp.1', names(clim))])
pairs( clim[ , grep('su.1', names(clim))])
# species effects --------------------------------------------------------------------- # 
dat <- merge( ARTR, clim , by = c('Treatment', 'year', 'Period'))
dat <- subset(dat, Period == 'Historical')

sw_f <- as.formula(  survives ~ logarea +
                       Group +
                       W.ARTR + 
                       W.HECO + 
                       W.POSE + 
                       W.PSSP + 
                       VWC.sp.1 + 
                       VWC.sp.0 +
                       VWC.sp.l +
                       VWC.su.0 + 
                       VWC.su.l + 
                       T.sp.1 + 
                       T.sp.0 +
                       T.sp.l + 
                       (logarea|year) )

clim_f <- as.formula( survives ~ logarea + 
                        Group +
                        W.ARTR + 
                        W.HECO + 
                        W.POSE + 
                        W.PSSP + 
                        P.f.w.sp.1 + 
                        P.f.w.sp.0 +  
                        P.f.w.sp.l + 
                        P.su.0 + 
                        P.su.l + 
                        T.sp.1 + 
                        T.sp.0 +
                        T.sp.l + 
                        (logarea|year) )


ARTR1 <- glmer(data = dat, sw_f, family = 'binomial')

summary(ARTR1)

ARTR2 <- lmer(data = dat, clim_f, family = 'binomial'  )

summary(ARTR2)

AIC( ARTR1, ARTR2)

# 
dat <- merge( HECO, clim , by = c('Treatment', 'year', 'Period'))
dat <- subset( dat, Period == 'Historical' & year > 1926)

HECO1 <- glmer(data = dat, sw_f, family = 'binomial' )

summary(HECO1)

HECO2 <- glmer(data = dat, clim_f , family = 'binomial')

summary(HECO2)

AIC( HECO1, HECO2)

#

dat <- merge( POSE, clim , by = c('Treatment', 'year', 'Period'))
dat <- subset( dat, Period == 'Historical' & year > 1926)

POSE1 <- glmer(data = dat, sw_f, family = 'binomial' )

summary(POSE1)

POSE2 <- glmer(data = dat, clim_f, family = 'binomial'  )

summary(POSE2)

AIC( POSE1, POSE2)

#

dat <- merge( PSSP, clim , by = c('Treatment', 'year', 'Period'))
dat <- subset( dat, Period == 'Historical' & year > 1926)

PSSP1 <- glmer(data = dat, sw_f, family = 'binomial' )

summary(PSSP1)

PSSP2 <- glmer(data = dat, sw_f, family = 'binomial' )

summary(PSSP2)

AIC( PSSP1, PSSP2)

