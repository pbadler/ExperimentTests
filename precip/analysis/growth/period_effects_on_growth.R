#
# test for period effects 
# 

library(lme4)

ARTR <- readRDS('data/temp_data/ARTR_growth.RDS')

ARTR1 <- lmer(data = ARTR, logarea.t1 ~ logarea.t0 + Group + Period + W.ARTR + Group + (logarea.t0|year)  )

summary(ARTR1)

# 

HECO <- readRDS('data/temp_data/HECO_growth.RDS')

HECO1 <- lmer(data = HECO, logarea.t1 ~ logarea.t0 + Group + Period + W.HECO + Group + (logarea.t0|year)  )

summary(HECO1)

#

POSE <- readRDS('data/temp_data/POSE_growth.RDS')

POSE1 <- lmer(data = POSE, logarea.t1 ~ logarea.t0 + Group + Period + W.POSE + Group + (logarea.t0|year)  )

summary(POSE1)

#

PSSP <- readRDS('data/temp_data/PSSP_growth.RDS')

PSSP1 <- lmer(data = PSSP, logarea.t1 ~ logarea.t0 + Group + Period + W.PSSP + Group + (logarea.t0|year)  )

summary(PSSP1)
