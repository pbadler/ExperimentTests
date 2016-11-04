model_string <- "
data {
  int<lower=0> N;
  vector[2] y[N];
}
parameters {
  vector[2] mu;
  cov_matrix[2] Sigma;
}
model {
  for (n in 1:N)
    y[n] ~ multi_normal(mu, Sigma);
}
generated quantities {
  real<lower=-1,upper=1> rho;
  rho <- Sigma[1,2] / sqrt(Sigma[1,1] * Sigma[2,2]);
}
"

library(MASS)
set.seed(1234)
N <- 1000
mu <- c(-2, 2)
Sigma <- matrix(c(10,3,3,2),nrow=2)
y <- mvrnorm(N, mu, Sigma)

library(rstan)
sm <- stan_model(model_code = model_string)
sf <- sampling(sm, data=c("N","y"))
