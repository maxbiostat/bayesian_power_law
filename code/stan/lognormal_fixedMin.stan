data{
  int<lower=0> N;
  real<lower=0> x_min;
  real<lower=x_min> x[N];
  real mean_mu;
  real<lower=0> sd_mu;
}
parameters{
  real mu;
  real<lower=0.0> sigma;
}
model{
 for(i in 1:N) target += lognormal_lpdf(x[i] | mu, sigma); 
   target += N* (log(2) - log(erfc( (log(x_min)-mu)/(sqrt(2)*sigma) )) ); // truncation correction
  /* priors */
  target += normal_lpdf(mu | mean_mu, sd_mu) ;
  target += gamma_lpdf(sigma | 1, 1);
  // target += normal_lpdf(sigma | 0, 1);
  // target += -2*log(sigma); // Jeffreys (joint) prior
}
