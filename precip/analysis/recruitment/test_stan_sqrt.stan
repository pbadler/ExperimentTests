// Single species model with climate: includes climate + intraspecific effects
data{
  int<lower=0> N; // observations
  int<lower=0> Nspp; // number of species 
  int<lower=0> spp; // focal species id
  matrix[N, Nspp] parents1; // parents in plot
  matrix[N, Nspp] parents2; // parents in group
  
}parameters{
  real<lower=0, upper=1> u;
}
transformed parameters{
  matrix[N, Nspp] trueP1;
  matrix[N, Nspp] trueP2;
  vector[N] p1; 
  vector[N] p2;
  vector[N] trueP1vec;
  vector[N] trueP2vec;
  matrix[N, Nspp] out; 
  
  p1 <- parents1[, Nspp];
  p2 <- parents2[, Nspp];

  trueP1 <- parents1 + parents2;

  for(n in 1:N){
    for( j in 1:Nspp) { 
      out[n, Nspp] <- sqrt(trueP1[n, j]);
    }
    //print(trueP1[n,])
  }
  
}
model{
  // Priors

  // Likelihood
  //Y ~ normal(10, 2);
}
