rm(list  = ls())
library(rstan)
# simulate climate, competition, year and group effects ------------------------------------- # 

test_dat <- readRDS('data/temp_data/survival_data_lists_for_stan.RDS')[['ARTR']]

sig_a <- 1
sig_b1 <- 0.2
b1_mu  <- 0.8

pars <- list( 
  bg = c(1.5, runif(5, -1, 1)),
  b2 = c(2, rep(0, 7)),
  a  = rnorm(test_dat$nyrs, 0, sd = sig_a),
  b1 = rnorm(test_dat$nyrs, b1_mu, sd = sig_b1),
  w  = c(1, 0, 0, 0))

simulate_survival <- function( pars , test_dat ){ 
  
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
  
  p <- exp(mu)/(1 + exp(mu)) # inverse logit 
  
  Y <- rbinom(N, size = 1,  p)
  rm(pars)
  detach(pars )   
  return(Y)
}

test_dat$Y <- simulate_survival(pars, test_dat)

load('data/temp_data/master_list.Rdata')

test_dat$tau_beta <- master_list$sd_vec[20, 'sd']
test_dat$tau_beta

myfit1 <- stan('analysis/survival/model_survival_1.stan', data = test_dat, chains = 4, cores = 4, iter = 1000)

traceplot(myfit1, 'w')

source('analysis/waic_fxns.R')

waic(myfit1)

ests1 <- summary(myfit1, c('b1_mu', 'w', 'b2', 'sig_a', 'sig_b1'))$summary[, 1]
ests1

mu <- summary( myfit1, 'mu')$summary[, 1]
muhat <- summary(myfit1, 'muhat')$summary[, 1]

plot( mu, test_dat$Y )
abline(0,1)

plot(muhat, test_dat$Yhold)
abline(0,1)

plot(test_dat$Xhold, muhat)



traceplot( myfit1, 'b2')
traceplot(myfit1, 'w')

# fit same with lmer 
library(lme4)

df <- data.frame( X = test_dat$X, Y = test_dat$Y)
df$C <- test_dat$C
df$W <- test_dat$W
df$Group <- factor(test_dat$gid)
df$yid <- test_dat$yid

f <- as.formula( Y ~ Group + X + (X|yid) + W + C)
m1 <- glmer(data = df, f , family = 'binomial')

summary(m1)

data.frame( summary(m1)$coefficients) 

Ceffects <- summary(myfit1, c('b2'))$summary[, 1]

cbind(fixef(m1)[grep( '^C' , names( fixef(m1)))], Ceffects, real = pars$b2)

# estimates table 

trueVals <- data.frame( type = 'true_val' , b1_mu = b1_mu, sig_b1 = sig_b1, sig_a = sig_a, 
                        Int = pars$bg[1], 
                        w_1 = pars$w[1], w_2 = pars$w[2], w_3 =pars$w[3], w_4 = pars$w[4], 
                        b2_1 = pars$b2[1], b2_2 = pars$b2[2], b2_3 = pars$b2[3], b2_4 = pars$b2[4], 
                        b2_5 = pars$b2[5], b2_6 = pars$b2[6], b2_7 = pars$b2[7], b2_8 = pars$b2[8]) %>% 
  gather(par , set_val , b1_mu:b2_7)

stanVals <- data.frame( summary(myfit1, c('bg[1]', 'b1_mu', 'sig_b1', 'sig_a', 'w', 'b2'))$summary[ , c(1,2)] )
stanVals$par <- row.names(stanVals) 

stanVals$par[ stanVals$par == "bg[1]" ] <- 'Int'
library(stringr)
stanVals$par <- str_replace_all( stanVals$par, pattern = '\\[', replacement = '_')
stanVals$par <- str_replace_all( stanVals$par, pattern = '\\]', replacement = '')

lmer_vals <- data.frame( summary(m1)$coefficients ) 
lmer_vals$par <- row.names(lmer_vals)

lmer_vals$par[ lmer_vals$par == 'X' ]  <- 'b1_mu'
lmer_vals$par[ grep('^WW',  lmer_vals$par) ]  <- paste0( 'w_', c(1:4))
lmer_vals

lmer_vals$par[ grep( '^C', lmer_vals$par ) ] <- paste0('b2_', c(1:8))

lmer_vals$par[lmer_vals$par == '(Intercept)'] <- 'Int'

merge( merge(trueVals, stanVals), lmer_vals, all.x = TRUE)

plot( data = merge( lmer_vals,  stanVals )[-1, ], Estimate ~ mean ) 
