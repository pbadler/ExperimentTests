ncores <- 4 
niter <- 2000 
nchains <- 4 
nthin <- 4

dl <- readRDS('output/stan_fits/ARTR_survival_model_data.RDS')

dl$X %>% head

X <- dl$X[, 1:5]
C1 <- dl$X[, 6]
C2 <- dl$X[, 7]

dl$X <- X
dl$C1 <- C2

dl$K <- ncol(dl$X)

mod <- rstan::stan_model('analysis/survival/logistic_gaussian2.stan') # load stan model 

fit <- rstan::sampling(mod, 
                       data = dl, 
                       chains = nchains, 
                       iter = niter, 
                       cores = ncores,
                       thin = nthin, 
                       pars = c('beta', 'center', 'scale', 'mu'), 
                       control = list(adapt_delta = 0.98))

center <- summary( fit, 'center')$summary 
scale <- summary( fit, 'scale')$summary

gauss <- function(x, center, scale){ 
  
  (1/(scale*sqrt(2*pi) ) )*exp( (-(1/2)*(( x - center )/scale)^2)) 
  
}
x <- seq(-2, 3, by = 0.1)
gauss(x, center = 0, scale = 2)
plot( x, gauss(x, center = 0, scale = 1) , type = 'l')

X <- dl$C1
S <- dl$S
S <- S[order(X)]
X <- X[order(X)]

mu <- inv_logit(gauss(X, center = center[,1], scale = scale[,1]) )
plot(X, S)
points(X, mu, type = 'l')
plot(X, mu, type = 'l')

