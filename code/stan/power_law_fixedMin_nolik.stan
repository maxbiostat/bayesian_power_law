functions{
  real power_law_lpdf(real x, real xm, real a) {
    real dens;
    // if(x < xm) return(negative_infinity());
     dens = log(a-1)-log(xm) -a *( log(x) - log(xm));
    return(dens);
  }
}
data{
  int<lower=0> N;
  real<lower=0> x_min;
  real<lower=x_min> x[N];
}
parameters{
  real<lower=1.0> alpha;
}
model{
 for(i in 1:N) target += power_law_lpdf(x[i] | x_min, alpha); 
  /* priors */
  target += -log(alpha - 1); // Jeffreys' prior
}
