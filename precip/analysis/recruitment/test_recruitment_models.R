rm(list  = ls())
library(rstan)

# simulate climate, competition, year and group effects ------------------------------------- # 
rm(list = ls())

test_dat <- readRDS('data/temp_data/modified_recruitment_data_lists_for_stan.RDS')[['ARTR']]

#pars <- rstan::extract(readRDS('output/stan_fits/predictions/PSSP_recruitment_1_19_4_predict.RDS'), c('bg', 'b2', 'sig_a', 'u', 'w', 'theta'))
sig_a <- 0.3

pars <- list(
  bg = c(1.5, runif(5, -1, 1)),
  b2 = c(1, rep(0, ncol(test_dat$C) - 1)),
  a  = rnorm(test_dat$nyrs, 0, sd = sig_a),
  u  = 0.8,
  w  = seq(0, 0, length.out = 4),
  theta = 1.1)

simulate_recruitment <- function( pars , test_dat ){ 
  
  nyrs <- test_dat$nyrs
  N <- test_dat$N
  C <- test_dat$C
  parents1 <- test_dat$parents1
  parents2 <- test_dat$parents2
  gm <- test_dat$gm
  yid <- test_dat$yid
  spp <- test_dat$spp
  Nspp <- test_dat$Nspp

  attach(pars )
  
  trueP1 <- parents1*u + parents2*(1-u)
  trueP2 <- sqrt(trueP1)

  climEff  <- C%*%b2
  gint     <- gm%*%bg
  coverEff <- trueP2%*%w
  
  mu <- coverEff
  lambda <- coverEff

  for(n in 1:N){
    mu[n] <- exp(gint[n] + a[yid[n]] + coverEff[n] + climEff[n]);
    lambda[n] <- trueP1[n, spp]*mu[n];  
  }
  Y <- rnbinom(N, mu = lambda, size = theta )
  detach(pars )   
  return(Y)
}

test_dat$Y <- simulate_recruitment( pars, test_dat )

myfit <- stan('analysis/recruitment/model_recruitment_1.stan', data = test_dat, iter = 2000 , cores = 4, seed = 1)

my_fit_summary <- summary(myfit, c('a', 'sig_a', 'b2'))$summary

plot( my_fit_summary[1:21, 1], pars$a)
abline(0,1)


myfit_ye <- stan('analysis/recruitment/model_recruitment_year_effects.stan', data = test_dat, iter = 2000 , cores = 4, seed = 1)
myfit_test <- stan('analysis/recruitment/model_recruitment_1_test.stan', data = test_dat, iter = 2000 , cores = 4, seed = 1)

r <- rstan::extract (myfit, 'lambda_pred')$lambda_pred
r_old <- rstan::extract( myfit, 'lambda')$lambda

new_df <- data.frame( obs =  test_dat$Yhold, predicted = colMeans(r) ) 
plot( data = new_df, obs ~ predicted)

old_df <- data.frame( obs = test_dat$Y, predicted = colMeans(r_old))
plot( data = old_df, obs ~ predicted)

rye <- rstan::extract (myfit_ye, 'lambda_pred')$lambda_pred
r_oldye <- rstan::extract( myfit_ye, 'lambda')$lambda

new_df_ye <- data.frame( obs =  test_dat$Yhold, predicted = colMeans(rye) ) 
plot( data = new_df_ye, obs ~ predicted)

old_df_ye <- data.frame( obs = test_dat$Y, predicted = colMeans(r_oldye))
plot( data = old_df_ye, obs ~ predicted)

rt <- rstan::extract (myfit_test, 'lambda_pred')$lambda_pred
r_oldt <- rstan::extract( myfit_test, 'lambda')$lambda

sig_a_test <- rstan::extract( myfit_test, 'sig_a')$sig_a
sig_a <- rstan::extract(myfit, 'sig_a')$sig_a
sig_a_ye <- rstan::extract(myfit_ye, 'sig_a')$sig_a

median(sig_a_ye)
median(sig_a_test)
median(sig_a)

new_dft <- data.frame( obs =  test_dat$Yhold, predicted = colMeans(rt) ) 

plot( data = new_dft , obs ~ predicted)
points( data = new_df, obs ~ predicted , col = 'red')

old_dft <- data.frame( obs = test_dat$Y, predicted = colMeans(r_oldt))
plot( data = old_dft, obs ~ predicted)
points( data = old_df, obs ~ predicted, col = 'red')







max(rye)
max(r)
max(r_oldye)
max(r_old)

new_df

mean( exp( rnorm(100, 0, 1)) )
mean( exp( rnorm(100, 0, 2)) )


var(colMeans(r_old))
var(colMeans(r))

hist( test_dat$C[,1])
hist( test_dat$Chold[,1])



data.frame( test_dat$yearhold, test_dat$Chold )

# 
# estimates <- summary(myfit, c('u', 'theta', 'sig_a', 'w', 'bg', 'a'))$summary[, 1]
# 
traceplot( myfit, 'b2')
traceplot( myfit, 'w')
# 
lambda <- summary(myfit, 'lambda')$summary[, 1]
plot( lambda, test_dat$Y ) 


# 
# 
# estimates[grep('b2', names(estimates))]
# plot ( pars$b2, estimates[grep('b2', names(estimates)) ])
# plot ( pars$bg, estimates[grep('bg', names(estimates))])
# plot(pars$a, estimates[grep('^a', names(estimates))])
# 
# 
# 
