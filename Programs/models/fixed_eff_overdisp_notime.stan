/* Simple linear regression */
data {
  int<lower=0> A; //number of age groups
  int<lower=0> R; //number of regions
  int<lower=0> N; //number of data points
  int<lower=0> y[A*R]; //vector of migrants
  vector[A*R] pop; //vector of log population
  matrix[N,(A-1)] age; //matrix of ages 
  matrix[N,(R-1)] region; //matrix of regions
}
parameters {
  vector[A-1] tau; //age-specific log rates
  vector[R-1] gamma; //region-specific intercept for log rates
  real alpha;// intercept
  vector[A*R] delta;
  real<lower=0> sigma; //sd of delta
}

transformed parameters{
  vector[N] log_rate_fixed; //vector of age & region specific log rate
  vector[N] log_rate; //vector of age & region specific log rate
      log_rate_fixed = alpha + age*tau + region*gamma ; 
      log_rate = alpha + age*tau + region*gamma + delta ; 
  }
model {
  //likelihood
      y ~ poisson(exp(log_rate + pop));
    
  //Constraint
  sum(delta) ~ normal(0, 0.01);

  //priors
  tau ~ normal(0,1);
  gamma ~ normal(0,1);
  alpha ~ normal(0,1);
  delta ~ normal(0, sigma);
  sigma ~ normal(0, 1);
  
  }
