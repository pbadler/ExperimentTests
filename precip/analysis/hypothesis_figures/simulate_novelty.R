#load libraries
library(rstan)

rm( list = ls() )

par(mfrow = c(1,1))

nsims = 100
#the explanatory variables
dat<-data.frame(x1=runif(nsims,-1,1), type = 'training')
dat$x2 <- dat$x1^2 # add quadratic effect if you want to simulate curve
dat <- dat[ order(dat$x1), ]

#the model
X<-model.matrix(~ x1 + x2,dat)

#the regression slopes
betas<-c(1, 2, 0) # an extra beta can be added to simulate a quadratic effect. Set to zero for linear model. 

#the standard deviation for the simulated data
sigma<-1

#the simulated training data
dat$y <- rnorm(nrow(X),X%*%betas,sigma)

# a matrix to get the predicted y values
newdat <- data.frame( x1 = runif( nsims, -1, 20))
newdat$x2 <- newdat$x1^2
newdat <- newdat[ order(newdat$x1), ]

new_X<-model.matrix(~ x1 + x2, newdat)

newdat$y2 <- rnorm( nrow(new_X), new_X%*%betas, sigma )

mydat <- list( new_X = new_X[, 1:2], 
               X = X[, 1:2], 
               K = ncol(X[, 1:2]), 
               y = dat$y, 
               y2 = newdat$y2, 
               N = length( dat$y ), 
               N2 = length(newdat$y2) )

plot( data = newdat, y2 ~ x1 , type = 'n') 
points(data = dat, y ~ x1, col = 'black'  )
points(data = newdat, y2 ~ x1, col = 'red'  )


# load the stan model from the e-mail 
myfit <- stan('analysis/hypothesis_figures/stan_lm.stan', data = mydat) ### stan model 

# stan summary 
dat$mu <- summary(myfit, 'linpred')$summary[,1]
newdat$mu <- summary(myfit, 'y_pred')$summary[, 1]
newdat$uci <- summary(myfit, 'y_pred')$summary[, 8]
newdat$lci <- summary(myfit, 'y_pred')$summary[, 4]

plot( data = newdat, y2 ~ x1 , type = 'n') 
points(data = dat, y ~ x1, col = 'black'  )
points(data = newdat, y2 ~ x1, col = 'red'  )

points( data = dat, mu ~ x1 , type = 'l')
points( data = newdat, mu ~ x1 , type = 'l', lty = 2)
points( data = newdat, uci ~ x1, type = 'l', col = 'blue', lty = 2)
points( data = newdat, lci ~ x1, type = 'l', col = 'blue', lty  = 2)

log_lik2 <- extract( myfit, 'log_lik2')$log_lik2
lpd <- log( colMeans(exp(log_lik2)))  # formula for log posterior predictive density

newdat$lpd <- lpd 
#
# X novelty defined as distance of new X from mean of X in training dat 
#
newdat$novelty <- abs(newdat$x2 - mean(X[,2])) 

#
# plot lpd and root squared error by novelty 
#
m1 = lm(data = newdat, lpd ~ novelty)
par(mfrow = c(1,2))
plot( data = newdat, lpd ~ novelty, col = 'red')
abline(m1, lwd = 2)

newdat$RSE <- sqrt( (newdat$y2 - newdat$mu)^2 ) # square root of squared prediction error
plot( data = newdat , RSE ~ novelty ) 
m2 = lm(data = newdat, RSE ~ novelty)
summary(m2)
abline(m2,  lwd= 2)

