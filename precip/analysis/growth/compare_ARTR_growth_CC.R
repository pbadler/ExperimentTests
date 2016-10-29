library(lme4)
library(dplyr)
library(tidyr)
rm(list = ls())

setwd('~/Documents/ExperimentTests/precip/')
ccmod <- readRDS('data/temp_data/cc_ARTR_growth_model.RDS')
mydat <- readRDS('data/temp_data/ARTR_scaled_growth_dataframe.RDS')
ccdat <- ccmod@frame

ccout <- lmer(logarea.t1~logarea.t0 + W + pptLag + ppt1 + TmeanSpr1 + 
           ppt2 + TmeanSpr2 + logarea.t0:pptLag + logarea.t0:ppt1 + 
           logarea.t0:TmeanSpr1 + logarea.t0:ppt2 + ppt1:TmeanSpr1 + 
           ppt2:TmeanSpr2 + logarea.t0:ppt1:TmeanSpr1 +
           (1|Group)+(logarea.t0|year),data=ccdat)      ##full model

summary(ccout)

# -------------my data ----------------------------------------------------- # 

ccdat_raw <- read.csv('data/temp_data/cc_ARTR_growth.csv') 
ccWdat_raw <- read.csv('data/temp_data/cc_ARTR_Ws_growth.csv')

mydat <- subset(mydat, Period == 'Historical')

ccdat_raw$year <- ccdat_raw$year + 1900
ccdat_raw <- cbind( ccdat_raw , ccWdat_raw ) 

ccdat_filtered <- merge( mydat[, c('year', 'trackID', 'quad', 'Group')], ccdat_raw, by = c('year', 'trackID', 'quad', 'Group' ))

ccdat_filtered <- ccdat_filtered %>% 
  arrange( year, Group, quad, trackID)

nrow(ccdat_filtered)
nrow(mydat)

W <- as.matrix( ccdat_filtered[ , grep('V', names(ccdat_filtered))] ) 

ccdat_filtered$W <- W

ccfiltered_out <- lmer(logarea.t1~logarea.t0 + W + pptLag + ppt1 + TmeanSpr1 + 
                         ppt2 + TmeanSpr2 + logarea.t0:pptLag + logarea.t0:ppt1 + 
                         logarea.t0:TmeanSpr1 + logarea.t0:ppt2 + ppt1:TmeanSpr1 + 
                         ppt2:TmeanSpr2 + logarea.t0:ppt1:TmeanSpr1 +
                         (1|Group)+(logarea.t0|year),data=ccdat_filtered)      ##full model

# rename my variables 

mydat$TmeanSpr1 <- mydat$T.sp.0
mydat$TmeanSpr2 <- mydat$T.sp.1
mydat$ppt1 <- mydat$P.f.w.sp.0
mydat$ppt2 <- mydat$P.f.w.sp.1
mydat$pptLag <- mydat$P.a.l

mydat <- 
  mydat %>% 
  arrange(year, Group, quad, trackID )

W <- as.matrix ( mydat[ , c('W.ARTR', 'W.HECO', 'W.POSE', 'W.PSSP')] ) 

mydat$W <- W
mydat_run <- mydat [ , names(ccdat)]

names(ccdat_filtered)
names(mydat_run)

myout <- lmer(logarea.t1~logarea.t0 + W + pptLag + ppt1 + TmeanSpr1 + 
                ppt2 + TmeanSpr2 + logarea.t0:pptLag + logarea.t0:ppt1 + 
                logarea.t0:TmeanSpr1 + logarea.t0:ppt2 + ppt1:TmeanSpr1 + 
                ppt2:TmeanSpr2 + logarea.t0:ppt1:TmeanSpr1 +
                (1|Group)+(logarea.t0|year),data=mydat_run)      ##full model

# estimates don't match !!!!!!!!!!!! 
summary( myout ) 
summary( ccfiltered_out)
summary(ccout)

# but variables are the same #####################
pairs( cbind( mydat_run$TmeanSpr2, ccdat_filtered$TmeanSpr2 ))
pairs( cbind( mydat_run$ppt1, ccdat_filtered$ppt1))
pairs( cbind( mydat_run$ppt2, ccdat_filtered$ppt2))
pairs( cbind( mydat_run$pptLag, ccdat_filtered$pptLag))


# rescale CC dataframe variables ################# 
ccdat_filtered2 <- ccdat_filtered
ccdat_filtered2[ ,c('ppt1', 'ppt2', 'pptLag', 'TmeanSpr1', 'TmeanSpr2')] <- scale( ccdat_filtered[ ,c('ppt1', 'ppt2', 'pptLag', 'TmeanSpr1', 'TmeanSpr2')] ) 

ccfiltered_out2 <- lmer(logarea.t1~logarea.t0 + W + pptLag + ppt1 + TmeanSpr1 + 
                          ppt2 + TmeanSpr2 + logarea.t0:pptLag + logarea.t0:ppt1 + 
                          logarea.t0:TmeanSpr1 + logarea.t0:ppt2 + ppt1:TmeanSpr1 + 
                          ppt2:TmeanSpr2 + logarea.t0:ppt1:TmeanSpr1 +
                          (1|Group)+(logarea.t0|year),data=ccdat_filtered2)      ##full model

summary(myout)
summary(ccfiltered_out2)

plot( fixef(myout), fixef(ccfiltered_out2)) ## basically the same 
plot( fixef(myout), fixef(ccfiltered_out))  ## 
points(fixef(myout)[1], fixef(ccfiltered_out)[1], col = 'red', cex= 3)
plot(fixef(myout)[-1], fixef(ccfiltered_out)[-1])
plot(fixef(myout), fixef(ccout))

summary( lm(logarea.t1~logarea.t0 + W + pptLag + ppt1 + TmeanSpr1 + 
              ppt2 + TmeanSpr2 + logarea.t0:pptLag + logarea.t0:ppt1 + 
              logarea.t0:TmeanSpr1 + logarea.t0:ppt2 + ppt1:TmeanSpr1 + 
              ppt2:TmeanSpr2 + logarea.t0:ppt1:TmeanSpr1, data = mydat_run) ) 

summary( lm(logarea.t1~logarea.t0 + W + pptLag + ppt1 + TmeanSpr1 + 
              ppt2 + TmeanSpr2 + logarea.t0:pptLag + logarea.t0:ppt1 + 
              logarea.t0:TmeanSpr1 + logarea.t0:ppt2 + ppt1:TmeanSpr1 + 
              ppt2:TmeanSpr2 + logarea.t0:ppt1:TmeanSpr1, data = ccdat_filtered2) ) 
