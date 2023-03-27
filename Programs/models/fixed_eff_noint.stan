/* Simple linear regression */
data {
  int<lower=0> A; //number of age groups
  int<lower=0> R; //number of regions
  int<lower=0> N; //number of data points
  int<lower=0> T; //number of years
  int<lower=0> y[A*R*T]; //vector of migrants
  vector[A*R*T] pop; //vector of log population
  int<lower=1> age[A*R*T]; //vector of ages 
  int<lower=1> region[A*R*T]; //vector of regions
  int<lower=1> time[A*R*T]; //vector of years
}
parameters {
  vector[A] tau; //age-specific log rates
  vector[R] gamma; //region-specific intercept for log rates
  vector[T] beta; //time-specific log rates
}

transformed parameters{
  vector[N] log_rate; //vector of age & region specific log rate
  for (i in 1:N){
      log_rate[i] = tau[age[i]] + gamma[region[i]] + beta[time[i]]; //age-specific rate + deviation
      }
    }
  
model {
  //likelihood
 
  y ~ poisson(exp(log_rate + pop));
  

  //priors
  tau ~ normal(0,1);
  gamma ~ normal(0,1);
  beta ~ normal(0,1);
  
  }
