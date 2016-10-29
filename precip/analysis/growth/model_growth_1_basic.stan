data{
  // training datalist, historical observations 
  int<lower=0> N;             // observations
  vector[N] Y;  // observation vector
  vector[N] X;                // size vector
  int<lower=0> G;             // groups
  matrix[N, G] gm;            // group dummy variable matrix 
  int<lower=0> nyrs;          // years
  int<lower=0> yid[N];        // year id
  int<lower=0> Wcovs;         // number of crowding effects 
  vector[N] W;                // crowding matrix
}
parameters{
  // for training data model  
  real<lower=0> sigma; 
  vector[G] bg;                     // varying group effects with first group as intercept 
  vector[nyrs] a;
  real b1_mu;
  vector[nyrs] b1;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real w;
}
transformed parameters{
// for training data model  
  vector[N] gint; 
  real mu[N];
  vector[N] crowdEff;

  // for training data model -----------------------------------
  gint <- gm*bg;
  crowdEff <- W*w;
  
  for(n in 1:N){
    mu[n] <- gint[n] + a[yid[n]] + b1[yid[n]]*X[n] + crowdEff[n];
  }
  
}
model{
   // for training data model 
  // Priors
  sigma ~ cauchy(0, 5);
  bg ~ normal(0,10);
  b1_mu ~ normal(0,10);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  a ~ normal(0,10);
  b1 ~ normal(b1_mu,10);
  w ~ normal(0,10);

  // Likelihood
  Y ~ normal(mu, sigma);
  
}
generated quantities {
  // hold out predictions 
  vector[N] log_lik;          // vector for computing log pointwise predictive density  

  for(n in 1:N){
      log_lik[n] <- normal_log(Y[n], mu[n], sigma);
  }
}

