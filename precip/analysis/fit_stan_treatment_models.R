#rm(list = ls())

library(rstan)

df <- expand.grid(species = c('ARTR', 'HECO', 'POSE', 'PSSP'), vital_rate = c('growth', 'recruitment', 'survival'))

i = 1
source('analysis/waic_fxns.R')
source('analysis/stan_data_functions.R')

spp <- df$species[i]
vr  <- df$vital_rate[i]

dat_file <- 'data/temp_data/ARTR_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

unique( cbind(dat$year, dat$yid) )

hold <- c(26:30)

# pars ------------------- 
formC <- as.formula(~-1)
small <- -1
formZ = as.formula(~ size) 
formE = as.formula(~ size)
formX = as.formula(~ size + small + Group + W)
# ---------------------- 
process_data <- function(dat, small = 0, formX, formC, formZ, formE, center = T, ... ){
  
  C <- model.matrix(formC, dat)
  dat$C <- scale(C)
  dat$W <- scale(dat$W)
  dat$Group <- factor(dat$gid)
  
  dat$Y <- scale(dat$logarea.t1)
  
  dat$size <- scale(dat$logarea.t0)
  dat$small <- factor(dat$logarea.t0 < small)
  dat$size_raw <- exp(dat$logarea.t0)

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

dl <- process_data(dat = dat, small = 0, formX = formX, formC = formC, formZ = formZ, formE = formE, center = T, historical = T, hold )

dl <- left_censor(dl, U = exp(-1))

gmod <- rstan::stan_model('analysis/growth/model_growth_censored.stan')

gfit <- rstan::sampling(gmod, 
                        data = dl, 
                        chains = 4, 
                        iter = 2000, 
                        cores = 4, 
                        pars = c('beta', 'eta', 'Y_hat'), 
                        control = list(adapt_delta = 0.9))

saveRDS(gfit, '~/Desktop/ARTR_growth_fit.RDS')
gfit <- readRDS('~/Desktop/ARTR_growth_fit.RDS')

Y <- dl$Y
sso <- shinystan::launch_shinystan(gfit)

length(Y)
test$summary %>% head

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


