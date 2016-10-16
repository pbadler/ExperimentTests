// Intra-specific competition model for survival: includes intraspecific effects only
// fit to entire dataset historical and modern observations 
data{
  int<lower=0> N2; // all observations
  int<lower=0,upper=1> Y2[N2]; // observation vector
  int<lower=0> Yrs2; // all years
  int<lower=0> yid2[N2]; // year id
  int<lower=0> gid2[N2]; // group id
  vector[N2] X2; // size vector
  int<lower=0> G; // groups
  int<lower=0> Wcovs; // number of crowding effects
  matrix[N2,Wcovs] W2; // crowding matrix
  int<lower=0> spp; // focal species 

}
parameters{
  real a_mu2;
  vector[Yrs2] a2;
  real b1_mu2;
  vector[Yrs2] b12;
  real gint2[G];
  vector[Wcovs] w2;
  real<lower=0> sig_a2;
  real<lower=0> sig_b12;
  real<lower=0> sig_G2;
  
}
transformed parameters{
  real mu2[N2];
  vector[N2] crowdEff2;

  crowdEff2 <- W*w2;

  for(n in 1:N2){
    mu2[n] <- inv_logit(a2[yid2[n]] + gint2[gid2[n]] + b12[yid2[n]]*X2[n] + crowdEff2[n]);
  }

}
model{
  a_mu2 ~ normal(0,10);
  w2 ~ normal(0,10);
  b1_mu2 ~ normal(0,10);
  sig_a2 ~ cauchy(0,2);
  sig_b12 ~ cauchy(0,2);
  sig_G2 ~ cauchy(0,2);
  for(g in 1:G)
    gint2[g] ~ normal(0, sig_G2);
  for(y in 1:Yrs2){
    a2[y] ~ normal(a_mu2, sig_a2);
    b12[y] ~ normal(b1_mu2, sig_b12);
  }

  // Likelihood
  Y2 ~ binomial(1,mu2);

}

