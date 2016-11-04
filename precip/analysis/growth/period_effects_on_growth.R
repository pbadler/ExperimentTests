#
# test for period effects 
# 

library(lme4)

ARTR <- readRDS('data/temp_data/ARTR_growth.RDS')

ARTR <- split(ARTR, ARTR$Period)

ARTR1 <- lmer(data = ARTR$Historical, logarea.t1 ~ logarea.t0 +W.ARTR+  Group + (logarea.t0|year)  )
ARTR2 <- lmer(data = ARTR$Historical, logarea.t1 ~ logarea.t0*W.ARTR+  Group + (logarea.t0|year)  )
AIC(ARTR1, ARTR2)
summary(ARTR1)

ARTR1p <- lmer(data = ARTR$Historical, logarea.t1 ~ poly( logarea.t0, 2) + W.ARTR + Group + (logarea.t0|year)  )

AIC(ARTR1, ARTR1p)
summary(ARTR1p)

ARTR <- do.call( rbind, ARTR)

ARTR$predicted1 <- predict(ARTR1p, newdata = ARTR, re.form = NA)
ARTR$deviation <- ARTR$predicted1 - ARTR$logarea.t0

par(mfrow = c(1,2)) 
plot(data = subset( ARTR, Period == 'Historical'), predicted1 ~ logarea.t1, main = 'historical')  
abline( 0, 1, col = 'red')

plot(data = subset( ARTR, Period == 'Modern'), predicted1 ~ logarea.t1, main = 'modern')  
abline( 0, 1, col = 'red')


aggregate( data = ARTR, deviation ~ Period , 'mean')

# --------------------------------------------------------------------- #

HECO <- readRDS('data/temp_data/HECO_growth.RDS')

HECO <- split(HECO, HECO$Period)

HECO1 <- lmer(data = HECO$Historical, logarea.t1 ~ logarea.t0 + W.HECO + Group + (logarea.t0|year)  )
HECO2 <- lmer(data = HECO$Historical, logarea.t1 ~ logarea.t0*W.HECO + Group + (logarea.t0|year))
AIC(HECO1, HECO2)
summary(HECO1)

HECO <- do.call( rbind, HECO)

HECO$predicted1 <- predict(HECO1, newdata = HECO, re.form = NA)
HECO$deviation <- HECO$predicted1 - HECO$logarea.t0

par(mfrow = c(1,2)) 
plot(data = subset( HECO, Period == 'Historical'), predicted1 ~ logarea.t1, main = 'historical')  
abline( 0, 1, col = 'red')

plot(data = subset( HECO, Period == 'Modern'), predicted1 ~ logarea.t1, main = 'modern')  
abline( 0, 1, col = 'red')


aggregate( data = HECO, deviation ~ Period , 'mean')

# --------------------------------------------------------------------- #

POSE <- readRDS('data/temp_data/POSE_growth.RDS')

POSE <- split(POSE, POSE$Period)

POSE1 <- lmer(data = POSE$Historical, logarea.t1 ~ logarea.t0+ W.POSE + Group + (logarea.t0|year)  )
POSE2 <- lmer(data = POSE$Historical, logarea.t1 ~ logarea.t0*W.POSE + Group + (logarea.t0|year)  )

POSE1p <- lmer(data = POSE$Historical, logarea.t1 ~ poly(logarea.t0,2) + W.POSE + Group + (logarea.t0|year)  )

AIC(POSE1, POSE2)
AIC(POSE1, POSE1p)
summary(POSE1p)

POSE <- do.call( rbind, POSE)

POSE$predicted1 <- predict(POSE1, newdata = POSE, re.form = NA)
POSE$deviation <- POSE$predicted1 - POSE$logarea.t0

par(mfrow = c(1,2)) 
plot(data = subset( POSE, Period == 'Historical'), predicted1 ~ logarea.t1, main = 'historical')  
abline( 0, 1, col = 'red')

plot(data = subset( POSE, Period == 'Modern'), predicted1 ~ logarea.t1, main = 'modern')  
abline( 0, 1, col = 'red')

aggregate( data = POSE, deviation ~ Period , 'mean')

# --------------------------------------------------------------------- #

PSSP <- readRDS('data/temp_data/PSSP_growth.RDS')

PSSP <- split(PSSP, PSSP$Period)

PSSP1 <- lmer(data = PSSP$Historical, logarea.t1 ~ logarea.t0 + W.PSSP + Group + (logarea.t0|year)  )
PSSP2 <- lmer(data = PSSP$Historical, logarea.t1 ~ logarea.t0*W.PSSP + Group + (logarea.t0|year)  )

AIC(PSSP1, PSSP2)
PSSP <- do.call( rbind, PSSP)

PSSP$predicted1 <- predict(PSSP1, newdata = PSSP, re.form = NA)
PSSP$deviation <- PSSP$predicted1 - PSSP$logarea.t0

par(mfrow = c(1,2)) 
plot(data = subset( PSSP, Period == 'Historical'), predicted1 ~ logarea.t1, main = 'historical')  
abline( 0, 1, col = 'red')

plot(data = subset( PSSP, Period == 'Modern'), predicted1 ~ logarea.t1, main = 'modern')  
abline( 0, 1, col = 'red')

aggregate( data = PSSP, deviation ~ Period , 'mean')

# --------------------------------------------------------------------- #
