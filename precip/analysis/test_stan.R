rm(list = ls())

load('data/temp_data/master_list.Rdata')

df <- readRDS('data/temp_data/growth_data_lists_for_stan.RDS')[['HECO']]
df$W3

library(rstan)
df$tau_beta <- master_list$sd_vec[1, 2]

# df$W <- df$W[, df$spp]
# df$W2 <- df$W2[, df$spp]
# df$Whold <- df$Whold[, df$spp]
# df$W3 <- df$W3[, df$spp]

testfit1 <- stan('analysis/growth/model_growth_3_predict.stan', chains = 1, data = df, cores = 1, iter = 100)



# init = rep(list(inits), 2) #

traceplot(testfit1, 'w')
traceplot(testfit1, 'bg')
traceplot(testfit1, 'sig_a')
traceplot(testfit1, 'theta')
traceplot(testfit1, 'u')

traceplot(testfit1, 'b2')

summary(testfit1, c('sig_a'))$summary
summary(testfit1, c('bg', 'a_raw', 'sig_a', 'w', 'sig_b1', 'b2'))$summary[, 1]


testfit1@model_pars

y_hat <- rstan::extract(testfit1, 'y_hat')$y_hat
q     <- rstan::extract(testfit1, 'q')$q
theta <- rstan::extract(testfit1, 'theta')$theta

q[1:10, 1]
mean(y_hat[, 1])
theta[1:10]


simulate_recruitment <- function(q, theta){ 
  out <- q
  out[]  <- NA
  for( i in 1:ncol(q)){ 
    out[, i] <-  mapply(function(mu, size) rnbinom(1, mu = mu, size = size ), mu = q[, i], size = theta)
  }
  return(out)
}

out <- simulate_recruitment(q, theta)

plot( colMeans(y_hat),  colMeans(out) ) 
max(out)
max(y_hat)

max( y_hat$y_hat ) 


summary(testfit1, c('b1_mu2', 'bg2', 'sig_a2', 'sig_b12', 'w2', 'sigma2'))$summary[, 1]

hist(df$W[df$W < 1])
max( df$W ) 
min( df$W)

# lmer model 

library(lme4)
spp_df <- readRDS('data/temp_data/POSE_scaled_growth_dataframe.RDS')
spp_df <- subset(spp_df, Period == 'Historical')
m1 <- lmer( data = spp_df, Y ~ gid + X + W.POSE + (X|yid)) 
summary(m1)
ranef(m1)

b1_mu = fixef(m1)['X']
b1 = ranef(m1)[[1]][, 2] + b1_mu

inits <- list( a = ranef(m1)[[1]][, 1], b1_mu = b1_mu, b1 = b1, bg = fixef(m1)[1:df$G], w = fixef(m1)['W.POSE'], sig_a = 1.5, sig_b1 = 0.25, tauSize = -0.5, tau = 4, sigma = 0.8)


plot(ranef(m1)[[1]][, 1], type ='l')

plot( summary(testfit1, c('a2'))$summary[, 6], type = 'l', col = 'red')
points( summary(testfit1, c('a'))$summary[, 6], type = 'l', col = 'blue')

plot(ranef(m1)[[1]][, 2], type ='l')
points( summary(testfit, c('b1'))$summary[, 6] - summary(testfit, 'b1_mu')$summary[, 1], type = 'l', col = 'red') 

test_preds <- rstan::extract(testfit, c('a', 'b1'))

test_df <- data.frame(rstan::extract(testfit1, c('y_hat')))
test_df2 <- data.frame(rstan::extract(testfit1, c('y_hat2')))

library(dplyr)
library(tidyr)
library(stringr)


y_hat_df <- test_df %>% gather( var, val , contains('y_hat')) %>% mutate( var_label = as.numeric(str_extract(pattern = '[0-9]+$', var)))
y_hat_df2 <- test_df2 %>% gather( var, val, contains('y_hat2')) %>% mutate(var_label = as.numeric(str_extract(pattern = '[0-9]+$', var)))

rdata <- readRDS('data/temp_data/PSSP_recruitment.RDS')
rdata_hold <- subset(rdata, Period == 'Modern')
nrow(rdata_hold)

rdata_hold$rowid <- 1:nrow(rdata_hold)

all_dat <- merge( y_hat_df, rdata_hold, by.x = 'var_label' , by.y = 'rowid')
all_dat2 <- merge( y_hat_df2, rdata_hold, by.x = 'var_label' , by.y = 'rowid')

all_dat <- all_dat %>% mutate( error = val - Y)
all_dat2 <- all_dat2 %>% mutate( error = val - Y)

ggplot (all_dat, aes( x = error)) + geom_histogram() + facet_wrap( ~ year, ncol = 4) + xlim(-20, 20)

ggplot (all_dat2, aes( x = error)) + geom_histogram() + facet_wrap( ~ year, ncol = 4) + xlim(-20, 20)


ggplot (a_df, aes( x = val)) + geom_histogram() + geom_vline(aes( xintercept = 0)) + facet_wrap( ~ var_label, ncol = 4)

