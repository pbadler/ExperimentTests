
f.ARTR <- as.formula(logarea.t1 ~  
                       W.ARTR + 
                       W.HECO + 
                       W.POSE + 
                       W.PSSP +
                       logarea.t0 + 
                       pptLag + logarea.t0:pptLag + 
                       ppt1 + logarea.t0:ppt1 + 
                       TmeanSpr1 + logarea.t0:TmeanSpr1 + 
                       ppt2 + logarea.t0:ppt2 + 
                       TmeanSpr2 +
                       ppt1:TmeanSpr1 +
                       ppt2:TmeanSpr2 + 
                       logarea.t0:ppt1:TmeanSpr1 + 
                       (logarea.t0|year) + 
                       (1|Group))

f.HECO <- as.formula(logarea.t1 ~ 
                       logarea.t0 + 
                       W.ARTR + 
                       W.HECO + 
                       W.POSE + 
                       W.PSSP +
                       pptLag + logarea.t0:pptLag + 
                       ppt1 + logarea.t0:ppt1 +     
                       TmeanSpr1 + logarea.t0:TmeanSpr1 + 
                       ppt2 + logarea.t0:ppt2 + 
                       TmeanSpr2 + logarea.t0:TmeanSpr2 + 
                       ppt1:TmeanSpr1 + 
                       ppt2:TmeanSpr2 + 
                      logarea.t0:ppt1:TmeanSpr1 + 
                  (logarea.t0|year) + 
                  (1|Group))

f.POSE <- as.formula(logarea.t1 ~ 
                       logarea.t0 + 
                       W.ARTR + 
                       W.HECO + 
                       W.POSE + 
                       W.PSSP +
                        pptLag + 
                        ppt1 +
                        TmeanSpr1 + 
                        ppt2  + 
                        TmeanSpr2 + 
                        logarea.t0:ppt1  + 
                        logarea.t0:TmeanSpr1  + 
                        logarea.t0:ppt2  + 
                        ppt1:TmeanSpr1  + 
                        ppt2:TmeanSpr2  + 
                        logarea.t0:ppt1:TmeanSpr1
                        (logarea.t0|year) + 
                        (1|Group))

f.PSSP <- as.formula(logarea.t1 ~ 
                       logarea.t0 + 
                       W.ARTR + 
                       W.HECO + 
                       W.POSE + 
                       W.PSSP +
                        pptLag + 
                        ppt1 + 
                        TmeanSpr1 + 
                        ppt2 + 
                        TmeanSpr2 + 
                        logarea.t0:pptLag + 
                        logarea.t0:ppt2 +
                        logarea.t0:TmeanSpr2 +
                        ppt2:TmeanSpr2 +
                        logarea.t0:ppt2:TmeanSpr2 
                        (logarea.t0|year) + 
                        (1|Group))

library(dplyr)
library(tidyr)

Cdat <- read.csv('~/driversdata/data/idaho/climateData/olderAdlerFiles/Climate.csv')
VWCdat <- readRDS('data/temp_data/all_clim_covs.RDS')
VWCdat <- subset( VWCdat, Period == 'Historical' & Treatment == 'Control')

Cdf <- merge( Cdat, VWCdat)

pairs(Cdf[, c('pptLag', 'ppt1', 'TmeanSpr1', 'ppt2', 'TmeanSpr2', 'VWC.sp.l', 'VWC.sp.0', 'VWC.sp.1', 'VWC.su.0', 'VWC.a.l')])
pairs(Cdf[, c('pptLag', 'VWC.a.l', 'P.a.l')])
pairs(Cdf[, c('ppt1', 'VWC.sp.0', 'P.f.w.sp.0')])

G <- readRDS('data/temp_data/HECO_growth.RDS')

df <- merge(G, Cdat, by = 'year')
df2 <- merge(G, VWCdat, by = 'year')

library(lme4)

myr <- lmer(data = df, 
                logarea.t1 ~ 
                  logarea.t0 + 
                  W.ARTR + 
                  W.HECO + 
                  W.POSE + 
                  W.PSSP + 
                  (logarea.t0|year) + 
                  (1|Group))

year_effects <- cbind( unique(df[, grep('^ppt', names(df))]), ranef(myr)$year)
VWC_year_effects <- cbind( unique(df2[, grep('^VWC', names(df2))]), ranef(myr)$year)
T_year_effects <- cbind( unique(df2[, grep('^T\\.', names(df2))]), ranef(myr)$year)

pairs(year_effects)
pairs( VWC_year_effects)
pairs(T_year_effects)

cor(year_effects)[, c('(Intercept)', 'logarea.t0')]
cor(T_year_effects)[, c('(Intercept)', 'logarea.t0')]
cor(VWC_year_effects)[ , c('(Intercept)', 'logarea.t0')]


cor.test(x = year_effects[, c('(Intercept)')], year_effects[, c('VWC.su.1')])
cor.test(x = year_effects[, c('(Intercept)')], year_effects[, c('VWC.su.0')])
cor.test(x = year_effects[, c('(Intercept)')], year_effects[, c('pptLag')])


m1 <- lmer(data = df, f.HECO)
m2 <- update(m1, . ~ .  - logarea.t0*ppt1*TmeanSpr1 - logarea.t0*ppt2*TmeanSpr2 + logarea.t0)
m3 <- update(m2, . ~ . - pptLag:logarea.t0)
m4 <- update(m3, . ~ . - pptLag)
summary(m1)
summary(m2)
summary(m3)
summary(m4)

AIC(m1, m2, m3, m4)





df_scaled <- df

df_scaled[, grep('^(Tmean)|^p', names(df_scaled)) ] <- scale(df_scaled[, grep('^(Tmean)|^p', names(df_scaled))])

m1s <- lmer(data = df_scaled, f.ARTR)

summary(m1)
summary(m1s)

plot(predict(m1), m1@frame$logarea.t1)
abline(0,1)

plot(predict(m1s), m1s@frame$logarea.t1)
abline(0,1)

library(glmmLasso)
df_scaled$year <- factor(df_scaled$year)
df_scaled

mlasso <- glmmLasso(logarea.t1 ~ logarea.t0 + ppt1, data = df_scaled, rnd = list(Group = ~1), lambda = 500)

plot(predict(mlasso), df_scaled$logarea.t1)
abline(0,1)

data("soccer")

soccer[,c(4,5,9:16)]<-scale(soccer[,c(4,5,9:16)],center=TRUE,scale=TRUE)
soccer<-data.frame(soccer)

## linear mixed model
lm1 <- glmmLasso(points ~ transfer.spendings + ave.unfair.score 
                 + ball.possession + tackles 
                 + ave.attend + sold.out, rnd = list(team=~1), 
                 lambda=500, data = soccer)

summary(lm1)

## similar linear model without random effects
lm1b <- glmmLasso(points ~ transfer.spendings + ave.unfair.score 
                  + ball.possession + tackles 
                  + ave.attend + sold.out, rnd = NULL, 
                  lambda=10, data = soccer)


plot(predict(lm1), soccer$points)
plot(predict(lm1b), df$logarea.t1)
length(predict(lm1b))
