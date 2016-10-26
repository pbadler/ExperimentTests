rm(list  = ls())
library(rstan)

# simulate climate, competition, year and group effects ------------------------------------- # 

test_dat <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')[['POSE']]

sig_a <- 1
sig_b1 <- 0.2
b1_mu  <- 0.8

pars <- list( 
  bg = c(1.5, runif(5, -1, 1)),
  b2 = c(1, rep(0, 13)),
  a  = rnorm(22, 0, sd = sig_a),
  b1 = rnorm(22, b1_mu, sd = sig_b1),
  w  = seq(0, 0, length.out = 4),
  sigma = 1.1)

simulate_recruitment <- function( pars , test_dat ){ 
  
  N <- test_dat$N
  C <- test_dat$C
  W <- test_dat$W
  X <- test_dat$X
  gm <- test_dat$gm
  yid <- test_dat$yid

  attach(pars )
  
  climEff  <- C%*%b2
  gint     <- gm%*%bg
  coverEff <- W%*%w
  
  mu <- coverEff
  for(n in 1:N){
    mu[n] <- gint[n] + a[yid[n]] + b1[yid[n]]*X[n] + coverEff[n] + climEff[n]
  }
  
  Y <- rnorm(N, mu, sigma)
  rm(pars)
  detach(pars )   
  return(Y)
}

test_dat$Y <- simulate_recruitment(pars, test_dat)
test_dat$W <- test_dat$W[, test_dat$spp]

hist( test_dat$X ) 
hist(test_dat$Y)

test_dat$tau_beta <- 10
myfit <- stan('analysis/growth/model_growth_2.stan', data = test_dat, iter = 1000, cores = 4 ) 

estimates <- summary(myfit, c('sig_a', 'w', 'bg', 'a', 'sigma', 'sig_b1', 'b2'))$summary[, 1]

plot ( pars$b2, estimates[grep('b2', names(estimates)) ])
plot ( pars$bg, estimates[grep('bg', names(estimates))])
plot(pars$a, estimates[grep('^a', names(estimates))])

estimates['b2[1]']
pars$b2
estimates['sigma']
pars$sigma 
estimates['sig_b1']
sig_b1
estimates['sig_a']
sig_a
