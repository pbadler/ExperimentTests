library(nls2)
setwd('~/Documents/ExperimentTests/precip/')
df <- readRDS('data/temp_data/POSE_growth.RDS')

df <- split(df, df$Period)

df$Historical$X <- exp( df$Historical$logarea.t0)
df$Historical$Y <- exp( df$Historical$logarea.t1)

plot(dat_asymp$logarea.t0, dat_asymp$logarea.t1)
abline( fit.lm, col = 'red')
abline(0, 1)

dat_asymp <- df$Historical

mean(dat_asymp$W.POSE)

fit.gomp2 <- nls(Y ~ a*exp(-b1*b2^X), start = list( a = 34, b1 = 2, b2 = 0.99 ), data = dat_asymp)
summary(fit.gomp2)

fit.gomp <- nls(Y ~ SSgompertz(X, Asym, b2, b3), trace = F, control = nls.control(maxiter=500), algorithm = "port", lower = c(Asym = 0, b2 = 0, b3 = 0), data = dat_asymp)
summary(fit.gomp)

fit.logis <- nls(Y ~ SSlogis(X, Asym, xmid, scal), trace = F, control = nls.control(maxiter=500), algorithm = "port", lower = c(Asym = 0, b2 = 0, b3 = 0), data = dat_asymp)
fit.lm <- lm( log(Y) ~ log( X )  , data = dat_asymp)

pred.gomp <- predict(fit.gomp, dat_asymp)
pred.logis <- predict(fit.logis, dat_asymp)

pred.lm   <- predict(fit.lm, dat_asymp)

plot(log( pred.gomp), log( dat_asymp$Y))
abline( 0, 1)

plot(log( pred.logis), log( dat_asymp$Y))
abline( 0, 1)

plot( pred.lm, log( dat_asymp$Y) ) 
abline( 0, 1)

m2 <- lm(log(Y) ~ poly( log(X), 2), data = dat_asymp)
plot( predict(m2), log(dat_asymp$Y))
abline( 0, 1)
summary(m2)

