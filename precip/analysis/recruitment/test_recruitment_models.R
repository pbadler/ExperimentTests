rm(list  = ls())
library(rstan)

# simulate recruitment model ----------------------------------------------------------------------------- # 
N <- 1000
X <- seq(1, 10, length.out = N)
beta0 <- 1 
beta1 <- 0.5
phi <- 0.5
mu <- exp( beta0 + beta1*X)
Y <- rnbinom(N, mu = mu, size = phi )
plot(X, Y)

test_dat <- list(Y = Y, N = N, X =X )

fit <- stan('analysis/recruitment/model_recruitment_simple.stan', data = test_dat, chains = 4)

summary(fit, c('beta1', 'beta0', 'phi'))$summary

traceplot(fit, 'phi' )

# with my data 

# simulate climate, competition, year and group effects ------------------------------------- # 
rm(list = ls())

test_dat <- readRDS('data/temp_data/recruitment_data_lists_for_stan.RDS')[['POSE']]

sig_a <- 1

pars <- list( 
  bg = c(1.5, runif(5, -1, 1)),
  b2 = c(0, 0, 0, 0, 0, 0, 0),
  a  = rnorm(22, 0, sd = sig_a),
  u  = 0.8,
  w  = seq(0, 0, length.out = 4),
  theta = 1.1)

simulate_recruitment <- function( pars , test_dat ){ 
  
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
  trueP2 <- trueP1
  for(n in 1:N)
    for( j in 1:Nspp)
      trueP2[n, j] <- sqrt(trueP1[n, j])
  
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

test_dat$Y <- simulate_recruitment(pars, test_dat)

test_dat$tau_beta <- 10
myfit <- stan('analysis/recruitment/model_recruitment_3.stan', data = test_dat, iter = 2000, cores = 4 ) 

estimates <- summary(myfit, c('u', 'theta', 'sig_a', 'w', 'bg', 'a'))$summary[, 1]

estimates[grep('b2', names(estimates))]
plot ( pars$b2, estimates[grep('b2', names(estimates)) ])
plot ( pars$bg, estimates[grep('bg', names(estimates))])
plot(pars$a, estimates[grep('^a', names(estimates))])

estimates['sig_a']
sig_a


