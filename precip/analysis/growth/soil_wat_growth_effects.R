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

ARTR <- readRDS('data/temp_data/ARTR_growth.RDS')
HECO <- readRDS('data/temp_data/HECO_growth.RDS')
POSE <- readRDS('data/temp_data/POSE_growth.RDS')
PSSP <- readRDS('data/temp_data/PSSP_growth.RDS')

# 
pairs( clim[, grep('sp.1', names(clim))])
pairs( clim[ , grep('su.1', names(clim))])
# species effects --------------------------------------------------------------------- # 
dat <- merge( ARTR, clim , by = c('Treatment', 'year', 'Period'))
dat <- subset(dat, Period == 'Historical')


sw_f <- as.formula(  logarea.t1 ~ logarea.t0 +
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
                       (logarea.t0|year) )

clim_f <- as.formula( logarea.t1 ~ logarea.t0 + 
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
                        (logarea.t0|year) )


sw_fp <- as.formula(  logarea.t1 ~ poly(logarea.t0, 2) +
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
                       (poly(logarea.t0, 2)|year) )

sw_fp1 <- as.formula(  logarea.t1 ~ poly(logarea.t0, 2) +
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
                        (logarea.t0|year) )

clim_fp <- as.formula( logarea.t1 ~ poly(logarea.t0, 2) + 
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
                        (poly(logarea.t0, 2)|year) )


fs <- list( sw_f, sw_fp, sw_fp1, clim_f, clim_fp) 

m <- lapply(fs , function(x) lmer(x, data = dat))

unlist (lapply(m, AIC)  ) 

# 
dat <- merge( HECO, clim , by = c('Treatment', 'year', 'Period'))
dat <- subset( dat, Period == 'Historical' & year > 1926)

m <- lapply(fs , function(x) lmer(x, data = dat))

unlist (lapply(m, AIC)  ) 

#

dat <- merge( POSE, clim , by = c('Treatment', 'year', 'Period'))
dat <- subset( dat, Period == 'Historical' & year > 1926)

m <- lapply(fs , function(x) lmer(x, data = dat))

unlist (lapply(mARTR, AIC)  ) 

#

dat <- merge( PSSP, clim , by = c('Treatment', 'year', 'Period'))
dat <- subset( dat, Period == 'Historical' & year > 1926)

m <- lapply(fs , function(x) lmer(x, data = dat))

unlist (lapply(mARTR, AIC)  ) 
