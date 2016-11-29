rm(list  = ls())
library(rstan)

# simulate climate, competition, year and group effects ------------------------------------- # 
rm(list = ls())

test_dat <- readRDS('data/temp_data/recruitment_data_lists_for_stan.RDS')[['ARTR']]

#pars <- rstan::extract(readRDS('output/stan_fits/predictions/PSSP_recruitment_1_19_4_predict.RDS'), c('bg', 'b2', 'sig_a', 'u', 'w', 'theta'))
sig_a <- 0.2

pars <- list(
  bg = c(1.5, runif(5, -1, 1)),
  b2 = c(1, rep(0, 7)),
  a  = rnorm(21, 0, sd = sig_a),
  u  = 0.8,
  w  = seq(0, 0, length.out = 4),
  theta = 1.1)
# test_dat$nyrs <- length(unique( test_dat$yid))

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
  
  a <- rnorm(nyrs, 0, sig_a)
  
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

test_dat$tau_beta <- 9

myfit <- stan('analysis/recruitment/model_recruitment_1.stan', data = test_dat, iter = 1000, cores = 4 )


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
