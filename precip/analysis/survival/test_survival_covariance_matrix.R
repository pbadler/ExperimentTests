rm(list = ls())

library(rstan)
library(tidyverse)
library(rstanarm)
library(lme4)
library(loo)

mod <- rstan::stan_model('analysis/survival/logistic.stan') # load survival model 

dat_file <- 'data/temp_data/ARTR_growth_survival_dataframe.RDS'
dat <- readRDS(dat_file)

dat <- dat %>% dplyr::select( survives, logarea.t0, W.ARTR, Period, year) %>% filter(Period == 'Historical')
S <- dat$survives
N <- length(S)
dat$logarea.t0 <- as.numeric(scale( dat$logarea.t0 ) )
dat$W.ARTR <- as.numeric(scale(dat$W.ARTR))

X <- model.matrix( ~ logarea.t0 + W.ARTR , data = dat)
g <- as.numeric(factor(dat$year))
G <- length(unique(g))
K <- ncol(X)
Z <- model.matrix( ~ logarea.t0 , data = dat)
J <- ncol(Z)

simple_dl <- list( S, N, X, g, G, K, Z, J )
names( simple_dl ) <- c('S', 'N', 'X', 'g', 'G', 'K', 'Z', 'J')

my_fit <- rstan::sampling(mod, data = simple_dl, cores = 4)

df <- data.frame(simple_dl)
df$X <- df$X.logarea.t0
df$W.ARTR <- df$X.W.ARTR

sg_fit <- stan_glmer('S ~ X + W.ARTR +  (X|g)' , data = df, cores = 4, family = 'binomial')

m1 <- glmer( 'S ~ X + W.ARTR + (X|g)', data = df , family = 'binomial')

Sigma_L <- summary(my_fit, 'Sigma_L')$summary
Sigma_L <- matrix( Sigma_L[, 1], 2,2, byrow = T)
SIGMA <- Sigma_L %*% t(Sigma_L)  
SIGMA

test_sg <- summary( sg_fit$stanfit, 'b')$summary
test_sg <- test_sg[1:52, 1]
sg_intercept <- test_sg[seq(1,52,by =2)]
sg_slope <- test_sg[1 + seq(1,52, by = 2)]
plot(sg_intercept, sg_slope)

test_mine <- summary( my_fit, 'u')$summary
test_mine <- test_mine[1:52, 1]
mine_intercept <- test_mine[seq(1,52, by = 2)]
mine_slope <- test_mine[1 + seq(1,52, by = 2)]
plot(mine_intercept, mine_slope)

plot(mine_intercept, sg_intercept)
abline(0,1)
plot(mine_slope, sg_slope)
abline(0,1)

SIGMA
as.matrix(bdiag(VarCorr(m1)))


