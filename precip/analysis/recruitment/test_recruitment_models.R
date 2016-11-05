rm(list  = ls())
library(rstan)

# simulate climate, competition, year and group effects ------------------------------------- # 
rm(list = ls())

test_dat <- readRDS('data/temp_data/recruitment_data_lists_for_stan.RDS')[['PSSP']]

sig_a <- 1

pars <- list( 
  bg = c(1.5, runif(5, -1, 1)),
  b2 = c(1, rep(0, 7)),
  a  = rnorm(21, 0, sd = sig_a),
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
  
  trueP2_scaled <- scale(trueP2)
  
  climEff  <- C%*%b2
  gint     <- gm%*%bg
  coverEff <- trueP2_scaled%*%w

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

test_dat$tau_beta <- 0.009

test_inits <- readRDS('data/temp_data/recruitment_init_vals.RDS')[['ARTR']]

inits <- rep( list ( test_inits), 4)

myfit <- stan('analysis/recruitment/model_recruitment_1.stan', data = test_dat, iter = 1000, cores = 4 ) 

estimates <- summary(myfit, c('u', 'theta', 'sig_a', 'w', 'bg', 'a'))$summary[, 1]

traceplot( myfit, 'w')
traceplot( myfit, 'lambda')

mu_pred <- summary(myfit, 'mu_pred')$summary


estimates[grep('b2', names(estimates))]
plot ( pars$b2, estimates[grep('b2', names(estimates)) ])
plot ( pars$bg, estimates[grep('bg', names(estimates))])
plot(pars$a, estimates[grep('^a', names(estimates))])



