data {
  int N; //the number of observations
  int N2; //the size of the new_X matrix
  int K; //the number of columns in the model matrix
  vector[N] y; //the response
  vector[N2] y2; // 
  matrix[N,K] X; //the model matrix
  matrix[N2,K] new_X; //the matrix for the predicted values
}
parameters {
  vector[K] beta; //the regression parameters
  real sigma; //the standard deviation
}
transformed parameters {
  vector[N] linpred;
  linpred <- X*beta;
}
model {  
  beta[1] ~ cauchy(0,10); //prior for the intercept following Gelman 2008
  for(i in 2:K)
  beta[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008

  y ~ normal(linpred,sigma);
}
generated quantities {
  vector[N2] y_pred;
  vector[N] log_lik;
  vector[N2] log_lik2;
  vector[N2] log_lik3; 
  y_pred <- new_X*beta; //the y values predicted by the model
  
  for( i in 1:N){
    log_lik[i] <- normal_log(y[i], linpred[i], sigma);
  }
  for( i in 1:N2){
    log_lik2[i] <- normal_log(y2[i], y_pred[i], sigma);
  }
  
  for ( i in 1:N2){ 
    log_lik3 [i ] <-normal_log(y_pred[i], y_pred[i], sigma);
  }
    
}
