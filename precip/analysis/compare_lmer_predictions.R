library(lme4)

spp <- 'HECO'

mydat <- read.csv(paste0 ('data/temp_data/', spp, '_growth_and_survival_cleaned_dataframe.csv'))

f0 <- 'logarea.t1 ~ logarea.t0 + W.ARTR + W.ARTR + W.POSE + W.PSSP + Group + (logarea.t0|year)'
f1 <- list( ARTR = 'logarea.t1 ~ logarea.t0 + W.ARTR + W.ARTR + W.POSE + W.PSSP + C.VWC.su.0 + C.T.sp.1 + Group + (logarea.t0|year)' , 
            HECO = 'logarea.t1 ~ logarea.t0 + W.ARTR + W.ARTR + W.POSE + W.PSSP + C.VWC.su.0 + C.VWC.su.1 + Group + (logarea.t0|year)' , 
            POSE = 'logarea.t1 ~ logarea.t0 + W.ARTR + W.ARTR + W.POSE + W.PSSP + C.T.su.l*logarea.t0 + Group + (logarea.t0|year)', 
            PSSP = 'logarea.t1 ~ logarea.t0 + W.ARTR + W.ARTR + W.POSE + W.PSSP + C.T.su.1 + C.T.f.1*lograea.t0 + Group + (logarea.t0|year)' )

fcc <- list( )
train <- subset(mydat, Period == 'Historical')
hold <- subset(mydat, Period == 'Modern')

m0 <- lmer(data= train, f0)
m1 <- lmer(data= train, f1[[spp]])

AIC(update( m0, REML = FALSE))
AIC(update( m1, REML = FALSE ))

hold$y_hat0 <- predict( m0, hold, re.form = NA)
hold$y_hat1 <- predict( m1, hold, re.form = NA)

y_hat_stan <- summary(readRDS(paste0('output/stan_fits/', spp, '_growth_climate_fit.RDS')), 'muhat')$summary[, 1]
standat <- readRDS('data/temp_data/modified_growth_data_lists_for_stan.RDS')[[spp]]
standat <- as.data.frame( standat[c('Xhold', 'Yhold', 'yearhold', 'yidhold', 'obs_idhold')] )
standat$obs_id <- standat$obs_idhold
standat$y_hat_stan <- y_hat_stan

hold <-  merge(hold, standat, by = 'obs_id') 

par(mfrow = c( 2, 2))
plot(hold$y_hat1, hold$y_hat_stan, main = paste('stan vs. lmer', spp )) # giving same predictions 
abline(0, 1) 

plot( hold$y_hat0, hold$logarea.t1, main = 'lmer year model')
abline(0,1, col = 'red')

plot( hold$y_hat1, hold$logarea.t1, main = 'lmer climate model')
abline(0,1, col = 'red')

plot( hold$y_hat_stan, hold$logarea.t1, main = 'stan climate model')
abline(0,1, col = 'red')

sum( resid(m0)^2 )/(sum( fitted(m0)^2))
sum( resid(m1)^2 )/(sum(fitted(m1)^2))

MSE0 <- mean( (hold$logarea.t1 - hold$y_hat0)^2)
MSE1 <- mean( (hold$logarea.t1 - hold$y_hat1)^2)
MSE_stan <- mean((hold$logarea.t1 - hold$y_hat_stan)^2 )
MSE0
MSE1
MSE_stan
summary(m0)

