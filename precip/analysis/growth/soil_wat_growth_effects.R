#
# test for period effects 
# 

library(lme4)

clim <- readRDS('data/temp_data/all_clim_covs.RDS')
clim[ , grep( '^[VPT]\\.', names(clim))] <- scale( clim[ , grep( '^[VPT]\\.', names(clim))])

ARTR <- readRDS('data/temp_data/ARTR_growth.RDS')
HECO <- readRDS('data/temp_data/HECO_growth.RDS')
POSE <- readRDS('data/temp_data/POSE_growth.RDS')
PSSP <- readRDS('data/temp_data/PSSP_growth.RDS')


# species effects --------------------------------------------------------------------- # 

dat <- merge( ARTR, clim , by = c('Treatment', 'year', 'Period'))
dat <- subset( dat, Period == 'Historical' & year > 1926)

ARTR1 <- lmer(data = dat, logarea.t1 ~ logarea.t0 + W.ARTR + Group + 
                VWC.sp.1_layer1 + 
                VWC.sp.1_layer2 + 
                VWC.sp.0_layer1 + 
                VWC.sp.0_layer2 + 
                T.sp.1 + 
                T.sp.0 + 
                (logarea.t0|year)  )

summary(ARTR1)

ARTR2 <- lmer(data = dat, logarea.t1 ~ logarea.t0 + W.ARTR + Group + 
                P.f.w.sp.1 + 
                P.f.w.sp.0 +  
                P.su.0 + 
                P.a.l + 
                T.sp.1 + 
                T.sp.0 + 
                (logarea.t0|year)  )

summary(ARTR2)

AIC( ARTR1, ARTR2)

# 

dat <- merge( HECO, clim , by = c('Treatment', 'year', 'Period'))
dat <- subset( dat, Period == 'Historical' & year > 1926)

HECO1 <- lmer(data = dat, logarea.t1 ~ logarea.t0 + W.HECO + Group + 
                VWC.sp.1_layer1 + 
                VWC.sp.1_layer2 + 
                VWC.sp.0_layer1 + 
                VWC.sp.0_layer2 + 
                VWC.su.0_layer1 + 
                VWC.su.0_layer2 + 
                VWC.a.l_layer1 + 
                VWC.a.l_layer2 + 
                T.sp.1 + 
                T.sp.0 + 
                (logarea.t0|year)  )

summary(HECO1)

HECO2 <- lmer(data = dat, logarea.t1 ~ logarea.t0 + W.HECO + Group + 
                P.f.w.sp.1 + 
                P.f.w.sp.0 +  
                P.su.0 + 
                P.a.l + 
                T.sp.1 + 
                T.sp.0 + 
                (logarea.t0|year)  )

summary(HECO2)

AIC( HECO1, HECO2)

#

dat <- merge( POSE, clim , by = c('Treatment', 'year', 'Period'))
dat <- subset( dat, Period == 'Historical' & year > 1926)

POSE1 <- lmer(data = dat, logarea.t1 ~ logarea.t0 + W.POSE + Group + 
                VWC.sp.1_layer1 + 
                VWC.sp.1_layer2 + 
                VWC.sp.0_layer1 + 
                VWC.sp.0_layer2 + 
                VWC.su.0_layer1 + 
                VWC.su.0_layer2 + 
                VWC.a.l_layer1 + 
                VWC.a.l_layer2 + 
                T.sp.1 + 
                T.sp.0 + 
                (logarea.t0|year)  )

summary(POSE1)

POSE2 <- lmer(data = dat, logarea.t1 ~ logarea.t0 + W.POSE + Group + 
                P.f.w.sp.1 + 
                P.f.w.sp.0 +  
                P.su.0 + 
                P.a.l + 
                T.sp.1 + 
                T.sp.0 + 
                (logarea.t0|year)  )

summary(POSE2)

AIC( POSE1, POSE2)

#

dat <- merge( PSSP, clim , by = c('Treatment', 'year', 'Period'))
dat <- subset( dat, Period == 'Historical' & year > 1926)

PSSP1 <- lmer(data = dat, logarea.t1 ~ logarea.t0 + W.PSSP + Group + 
                VWC.sp.1_layer1 + 
                VWC.sp.1_layer2 + 
                VWC.sp.0_layer1 + 
                VWC.sp.0_layer2 + 
                VWC.su.0_layer1 + 
                VWC.su.0_layer2 + 
                VWC.a.l_layer1 + 
                VWC.a.l_layer2 + 
                T.sp.1 + 
                T.sp.0 + 
                (logarea.t0|year)  )

summary(PSSP1)

PSSP2 <- lmer(data = dat, logarea.t1 ~ logarea.t0 + W.PSSP + Group + 
                P.f.w.sp.1 + 
                P.f.w.sp.0 +  
                P.su.0 + 
                P.a.l + 
                T.sp.1 + 
                T.sp.0 + 
                (logarea.t0|year)  )

summary(PSSP2)

AIC( PSSP1, PSSP2)
