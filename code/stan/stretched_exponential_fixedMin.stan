functions{
  real stretched_exponential_lpdf(real x, real beta, real lambda, real xmin){
    real lc = log(beta) + log(lambda) + lambda * xmin^beta;
    real ldens = (beta-1) * log(x) -lambda * x^beta;
    return(lc + ldens);
  }
}
data{
  int<lower=0> N;
  real<lower=0> x_min;
  real<lower=x_min> x[N];
  real<lower=0> a1;
  real<lower=0> b1;
  real<lower=0> a2;
  real<lower=0> b2;
}
parameters{
  real<lower=0> beta;
  real<lower=0> lambda;
}
model{
 for(i in 1:N) target += stretched_exponential_lpdf(x[i] | beta, lambda, x_min); 
  /* priors */
  target += gamma_lpdf(beta | a1, b1) ;
  target += gamma_lpdf(lambda | a2, b2);
}
