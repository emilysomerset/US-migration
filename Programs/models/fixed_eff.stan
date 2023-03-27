/* Simple linear regression */
data {
  int<lower=0> A; //number of age groups
  int<lower=0> R; //number of regions
  int<lower=0> N; //number of data points
  int<lower=0> T; //number of years
  int<lower=0> y[A*R*T]; //vector of migrants
  vector[A*R*T] pop; //vector of log population
  matrix[N,(A-1)] age; //matrix of ages 
  matrix[N,(R-1)] region; //matrix of regions
  matrix[N,(T-1)] time; //matrix of years
}
parameters {
  vector[A-1] tau; //age-specific log rates
  vector[R-1] gamma; //region-specific intercept for log rates
  vector[T-1] beta; //time-specific log rates
  real alpha;// intercept
}

transformed parameters{
  vector[N] log_rate; //vector of age & region specific log rate
      log_rate = alpha + age*tau + region*gamma + time*beta; 
  }
model {
  //likelihood
      y ~ poisson(exp(log_rate + pop));

  //priors
  tau ~ normal(0,1);
  gamma ~ normal(0,1);
  alpha ~ normal(0,1);
  beta ~ normal(0,1);
  
  }
