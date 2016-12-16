
setwd('~/Documents/ExperimentTests/precip/')
myclim <- readRDS('data/temp_data/all_clim_covs.RDS')
mydat <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')
ccdat <- read.csv('~/Desktop/Data & Scripts/wrapper-ID/climateData/Climate.csv')

head(ccdat)

df <- merge( myclim , ccdat, by = 'year')

plot(df$ppt2, df$P.f.w.sp.1)
abline(0,1)

plot(df$ppt1, df$P.f.w.sp.0)
abline(0,1)

plot(df$TmeanSpr1, df$T.sp.0)
abline(0,1)

myClim <- unique( data.frame( year = mydat$ARTR$year, mydat$ARTR$C ) )

head( myClim)
head(ccdat)

all_dat <- merge( myClim, ccdat, by = 'year')

head(all_dat)

plot(scale(all_dat$TmeanSpr2), all_dat$T.sp.1)
abline(0,1)
plot(scale(all_dat$ppt2), all_dat$VWC.sp.1)
abline(0,1)

all_dat <- read.csv('data/temp_data/ARTR_growth_and_survival_cleaned_dataframe.csv')
all_dat <- subset(all_dat, Period == 'Historical')

all_dat <- merge(all_dat, ccdat, by = 'year')

cdat <- unique( test[ , c('year', 'ppt1', 'ppt2', 'C.VWC.sp.0', 'C.VWC.sp.1', 'TmeanSpr1', 'TmeanSpr2', 'C.T.sp.0', 'C.T.sp.1') ] )
head(cdat)
plot(cdat$ppt1, cdat$C.VWC.sp.0)
plot(cdat$ppt2, cdat$C.VWC.sp.1)
plot(cdat$TmeanSpr1, cdat$C.T.sp.0)
plot(cdat$TmeanSpr2, cdat$C.T.sp.1)

m <- lmer(data = all_dat, 'logarea.t1 ~ X + W.ARTR + W.POSE + W.PSSP + W.HECO + C.VWC.su.0 + C.T.sp.1 + (1|Group) + (X|yid)')
summary(m)

mnull <- lmer(data = all_dat, 'logarea.t1 ~ X + W.ARTR + W.POSE + W.PSSP + W.HECO + (1|Group) + (X|yid)')
AIC(update(m, REML = F))
AIC(update(mnull, REML = F))

ye <- read.csv('output/ARTR_growth_year_effects_table.csv')

head(ye)

ye <- ye[1:21, 2]
cor(unique( mydat$C.T.sp.1)[1:21], ye)
cor(unique(mydat$C.T.sp.0)[1:21], ye)

ye
unique( mydat$C.T.sp.1mydat$C.T.su.0 )


ccmodel <- lmer(data = all_dat, 'logarea.t1 ~ X + W.ARTR + W.HECO + W.POSE + W.PSSP + (1|Group) + (X|yid) + 
        pptLag + ppt1 + TmeanSpr1 + 
         ppt2 + TmeanSpr2 + X:pptLag + X:ppt1 + 
         X:TmeanSpr1 + X:ppt2 + ppt1:TmeanSpr1 + 
         ppt2:TmeanSpr2 + X:ppt1:TmeanSpr1')



ccmodel2 <- lmer(data = all_dat, 'logarea.t1 ~ X + W.ARTR + W.HECO + W.POSE + W.PSSP + (1|Group) + (X|yid) + 
                  C.VWC.sp.l + C.VWC.sp.0 + 
                   C.T.sp.0 + C.VWC.sp.1 + C.T.sp.1 + X:C.VWC.sp.l + X:C.VWC.sp.0 + 
                   X:C.T.sp.0 + X:C.VWC.sp.1 + ppt1:C.T.sp.0 + 
                   ppt2:C.T.sp.1 + X:C.VWC.sp.0:C.T.sp.0')



summary(ccmodel)
summary(ccmodel2)
AIC(update(ccmodel, REML = F))
AIC( update(mnull, REML = F))


mydat <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')
ccClim <- merge(data.frame( year = mydat$ARTR$year),  ccdat, by = 'year', all.x = T)

as.matrix(ccClim)
test <- cbind(ccClim, mydat$ARTR$C)
par(mfrow = c(1,1))
plot(test$pptLag, test$VWC.sp.l)
plot(test$ppt1, test$VWC.sp.0)
plot(test$ppt2, test$VWC.sp.1)




mydat$ARTR$C <- scale( as.matrix(ccClim[, -1]))
mydat$ARTR$Covs <- ncol(mydat$ARTR$C)


m1 <- rstan::stan('analysis/growth/model_growth_test.stan', data = mydat$ARTR, chains = 4, cores =4)
mydat$ARTR$C <- mydat$ARTR$C[, -c(1,3)]
mydat$ARTR$Covs <- ncol(mydat$ARTR$C)
m1 <- rstan::stan(fit = m1, data = mydat$ARTR, chains = 4, cores = 4)

m2 <- rstan::stan('analysis/growth/model_growth_year_effects.stan', data = mydat$ARTR, chains = 4, cores = 4)

old_waic <- waic( readRDS('output/stan_fits/ARTR_growth_climate_fit.RDS') )


waic(m1)
waic(m2)
old_waic

traceplot(m1, 'b2')


mt_new <- rstan::stan('analysis/growth/model_growth_treatment_effects2.stan', data = mydat$ARTR, chains = 4, cores = 4)
mt_old <- rstan::stan('analysis/growth/model_growth_treatment_effects.stan', data = mydat$ARTR, chains = 4, cores = 4)

par(mfrow = c(1,2))
traceplot(mt_new, 'bt')
traceplot(mt_old, 'bt')

extract(mt_new, 'bt')$bt

waic(mt_new, 'log_lik2')

mydat$ARTR$tm2 <- mydat$ARTR$tm2[, -c(3:4)]
mydat$ARTR$nT <- 2
mydat$ARTR$tm3 <- mydat$ARTR$tm3[, -c(3:4)]
mt_new.2 <- rstan::stan('analysis/growth/model_growth_treatment_effects2.stan', data = mydat$ARTR, chains = 4, cores = 4)

waic(mt_new.2, 'log_lik2')
