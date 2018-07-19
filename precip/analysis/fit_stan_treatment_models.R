rm(list = ls())

library(rstan)

df <- expand.grid(species = c('ARTR', 'HECO', 'POSE', 'PSSP'), vital_rate = c('growth', 'recruitment', 'survival'))

i = 1
source('analysis/waic_fxns.R')
source('analysis/stan_data_functions.R')

spp <- df$species[i]
vr  <- df$vital_rate[i]

dat_file <- 'data/temp_data/POSE_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

process_data <- function(dat, formX, formC, formZ, formE, center = T, ... ){
  
  C <- model.matrix(formC, dat)
  dat$C <- scale(C)
  dat$W <- scale(dat$W)
  dat$Group <- factor(dat$gid)
  
  dat$X <- model.matrix(formX, data = dat)
  dat$Z <- model.matrix(formZ, data = dat)
  dat$E <- model.matrix(formE, data = dat) 
  
  dat$g <- factor(dat$yid)
  
  dat_4_cover <- dat ### Need to preserve dataframe with NA's (dead plants) for predicting cover 
  dat_4_cover <- split_df(dat_4_cover, hold = 0)
  dl_4_cover <- make_dl(dat_4_cover)
  dl_4_cover <- dl_4_cover[-grep('hold', names(dl_4_cover))]
  
  dat <- dat[complete.cases(dat), ]
  dat <- split_df(dat, hold )
  dl <- make_dl(dat)
  
  names(dl_4_cover) <- paste0( 'cover_', names(dl_4_cover))
  
  return( c(dl, dl_4_cover))
}


# pars ------------------- 
hold <- c(26:30)
formC <- as.formula(~-1)
small <- -1
formZ = as.formula(~ size) 
formE = as.formula(~ size)
formX = as.formula(~ size + W + GroupP2 + Treatment2)
# ---------------------- 

dat$size <- scale( dat$logarea.t0 )
dat$Y    <- scale( dat$logarea.t1 )
hist(dat$Y)
dat$small <- as.numeric(dat$size < -1)
dat$GroupP2 <- as.numeric( dat$Group == 'P2')

dl <- process_data(dat = dat, formX = formX, formC = formC, formZ = formZ, formE = formE, center = T, historical = T, hold )

freq <- data.frame( table(dl$X[,2]) )
freq[ order(freq$Var1), ][1:15, ] # look at the smallest plants to decide cutoff

dl <- left_censor(dl, U = -1.11)

gmod <- rstan::stan_model('analysis/growth/model_growth_censored.stan')

gfit1 <- rstan::sampling(gmod, 
                        data = dl, 
                        chains = 4, 
                        iter = 2000, 
                        cores = 4, 
                        pars = c('beta', 'eta', 'Y_hat', 'u', 'log_lik'))

saveRDS(gfit1, '~/Desktop/POSE_growth_fit1.RDS')

### posterior predictive check
Y_hat <- data.frame( summary( gfit1, 'Y_hat')$summary )
Y_hat$obs <- dl$Y
Y_hat$size <- dl$X[,2]
library(tidyverse)

mean( Y_hat$obs < Y_hat$X75. & Y_hat$obs > Y_hat$X25.  )
mean( Y_hat$obs < Y_hat$X97.5. & Y_hat$obs > Y_hat$X2.5.  )

mean( Y_hat$obs < Y_hat$X25. )

Y_hat[ sample( 1:nrow(Y_hat), 50 ), ]  %>% 
  ggplot( aes(x = obs, y = mean, ymin = X25., ymax = X75.)) + 
  geom_point( aes( y = obs), col = 'red') + 
  geom_errorbar(col = 'blue')

ll <- loo::extract_log_lik(gfit1)
loo1 <- loo::loo(ll)
loo1

Y <- dl$Y
colnames(dl$X)
sso <- shinystan::launch_shinystan(gfit1)

# Try a different model 
# pars ------------------- 
hold <- c(26:30)
formC <- as.formula(~-1)
small <- -1
formZ = as.formula(~ size) 
formE = as.formula(~ size)
formX = as.formula(~ size + small + W + GroupP2 + Treatment2)
# ---------------------- 

dl <- process_data(dat = dat, formX = formX, formC = formC, formZ = formZ, formE = formE, center = T, historical = T, hold )
dl <- left_censor(dl, U = -3)

gfit2 <- rstan::sampling(gmod, 
                        data = dl, 
                        chains = 4, 
                        iter = 2000, 
                        cores = 4, 
                        pars = c('beta', 'eta', 'Y_hat', 'u', 'log_lik'))

saveRDS(gfit2, '~/Desktop/ARTR_growth_fit2.RDS')

### posterior predictive check
Y_hat <- data.frame( summary( gfit2, 'Y_hat')$summary )
Y_hat$obs <- dl$Y
Y_hat$size <- dl$X[,2]

mean( Y_hat$obs < Y_hat$X75. & Y_hat$obs > Y_hat$X25.  )
mean( Y_hat$obs < Y_hat$X97.5. & Y_hat$obs > Y_hat$X2.5.  )

mean( Y_hat$obs < Y_hat$X25. )

Y_hat[ sample( 1:nrow(Y_hat), 50 ), ]  %>% 
  ggplot( aes(x = obs, y = mean, ymin = X25., ymax = X75.)) + 
  geom_point( aes( y = obs), col = 'red') + 
  geom_errorbar(col = 'blue')

ll <- loo::extract_log_lik(gfit2)
loo2 <- loo::loo(ll)

loo1
loo2


year_effects <- data.frame( summary(gfit1, 'u')$summary )

intercepts <- year_effects[ grep( ',1', row.names(year_effects)), ]
slopes <- year_effects[ grep( ',2', row.names(year_effects)), ]

plot( unique(dat$year)[1:25], intercepts$mean, type = 'l')
plot( unique(dat$year)[1:25], slopes$mean, type = 'l')

library(lme4)
m1 <- lmer(Y ~ size + W + GroupP2 + Treatment2 + (size|yid), data = dat[!dat$yid %in% hold, ])
m1_stan <- rstanarm::stan_glmer(Y ~ size + W + GroupP2 + Treatment2 + (size|yid), data = dat[!dat$yid %in% hold, ], cores = 4)

loo_stanarm1 <- loo::loo(m1_stan)
loo_stanarm1
loo1


plot( unique(dat$year)[1:25], ranef(m1)$yid[, 1], type = 'l')
plot( unique(dat$year)[1:25], ranef(m1)$yid[, 2], type = 'l')

plot( unique(dat$year)[1:25], ranef(m1_stan)$yid[, 1], type = 'l')
plot( unique(dat$year)[1:25], ranef(m1_stan)$yid[, 2], type = 'l')

plot( intercepts$mean, ranef(m1)$yid$`(Intercept)` )
plot( slopes$mean, ranef(m1)$yid$`size`)

plot( m1_stan$coefficients)

stanarm_year_effects <- data.frame( rstan::summary(m1_stan$stanfit, 'b')$summary )

stanarm_I <- stanarm_year_effects[ grep('Intercept', row.names(stanarm_year_effects)) , ] 
stanarm_slope <- stanarm_year_effects[ grep('size', row.names(stanarm_year_effects)) , ] 
plot(stanarm_I$mean[-26], ranef(m1_stan)$yid[,1])
plot(stanarm_slope$mean[-26], ranef(m1_stan)$yid[,2])




plot(Y, test$summary[,1])
test <- (summary(gfit, 'Y_hat'))

nrow(test$summary)

for(i in 1:nrow(df)){ 
  
  spp <- df$species[i]
  vr  <- df$vital_rate[i]
  
  dat <- readRDS(paste0('data/temp_data/', vr, '_data_lists_for_stan.RDS'))[[spp]] # fit full treatment by size effects
  
  # drop all size * treatment effects from design matrix
  dat$tm2 <- dat$tm2[, -c(3:4)] # remove size by treatment effects 
  dat$tm3 <- dat$tm3[, -c(3:4)] # remove size by treatment effects 
  dat$nT <- ncol(dat$tm2)
  browser()
  
  myfit <- stan(paste0('analysis/', vr, '/model_', vr, '_treatment_effects2.stan'), data = dat, cores = 4, iter = 2000, thin = 4, seed = 1)
  
  # if(vr != 'recruitment' ) { ### Deal with size by treatment effects 
  #   
  #   waic1 <- waic(myfit, 'log_lik2')
  #   
  #   dat2 <- dat
  #   dat2$tm2 <- dat2$tm2[, -c(3:4)] # remove size by treatment effects 
  #   dat2$tm3 <- dat2$tm3[, -c(3:4)] # remove size by treatment effects 
  #   dat2$nT <- ncol(dat2$tm2)
  #   
  #   myfit2 <- stan( fit = myfit, data = dat2, cores = 4, iter = 2000, thin = 4, seed = 1)
  #   
  #   waic2 <- waic(myfit2, 'log_lik2')
  #   
  #   # if removing the size by treatment parameters improves fit then take them out  -------------------------- # 
  #   if ( waic1$waic < waic2$waic ){ 
  #     myfit <- myfit 
  #     dat <- dat 
  #   }else if( waic1$waic > waic2$waic ) { 
  #     myfit <- myfit2  
  #     dat <- dat2 
  #   }
  # 
  # }
  
  
  # ----------------------------------------------------------------------------------------------------------#   
  
  ss <-  get_sampler_params(myfit) 
  
  dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
  
  
  if ( dv > 0 ) { 
    # try again if divergence 
    inits <- apply(myfit, 2, relist, skeleton = rstan:::create_skeleton(myfit@model_pars, myfit@par_dims)) # find initial values from end of chain 
    myfit <- stan( fit = myfit, init = inits, data = dat, cores = 4, iter = 4000, thin = 8, seed = 1)
  }
  
  ss <-  get_sampler_params(myfit) 
  
  dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
  
  if ( dv > 0 ){ 
    # try again if divergent 
    inits <- apply(myfit, 2, relist, skeleton = rstan:::create_skeleton(myfit@model_pars, myfit@par_dims))
    myfit <- stan( fit = myfit, init = inits, data = dat, cores = 4, iter = 2000, thin = 4, 
                  control = list(adapt_delta = 0.85, stepsize = 0.8, max_treedepth = 20), seed = 1 )
  } 
  
  ss <-  get_sampler_params(myfit) 
  
  dv <- sum(   unlist( lapply( ss, function(x) sum( x[ (1 + ceiling(0.5*nrow(x))):nrow(x), 'divergent__']) )))
  
  if ( dv > 0 ){ 
    # try again if divergent 
    inits <- apply(myfit, 2, relist, skeleton = rstan:::create_skeleton(myfit@model_pars, myfit@par_dims)) 
    myfit <- stan( fit = myfit, init = inits, data = dat, cores = 4, iter = 4000, thin = 8, 
                   control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 20), seed = 1 )
  } 
  
  saveRDS(myfit, paste0('output/stan_fits/', spp, '_', vr, '_treatment_fit.RDS'))
  
  rm(myfit, dat, inits)

}


