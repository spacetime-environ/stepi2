data {
  int<lower=0> N;
  real<lower=0> E[N]; // need to indicate that variable is strictly positive
  int<lower=0> Y[N];
}

parameters {
  real<lower=0> theta[N];
  real<lower=0> a;
  real beta0;
}

transformed parameters{
  real <lower=0> mu[N];
  
  for(i in 1:N){
    mu[i]=E[i]*theta[i]*exp(beta0);
  }
  
}

model {
  // likelihood function and prior for theta
  for(i in 1:N){
    Y[i] ~ poisson(mu[i]);
    theta[i]~gamma(a,a);
  }
  a~gamma(1, 1);
  beta0~normal(0,10);
}

generated quantities {
  vector [N] log_lik;
  int<lower=0> yfit [N];
  
  //computing the log_likelihood for each value of the mean mu and the fitted values
  for(i in 1:N){
    log_lik[i]=poisson_lpmf(Y[i] |mu[i]);
    yfit[i]=poisson_rng(mu[i]);
  }
  
}
