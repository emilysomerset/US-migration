/* RW on overdispersion */
data {
  int<lower=0> A; //number of age groups
  int<lower=0> R; //number of regions
  int<lower=0> N; //number of data points
  int<lower=0> y[A*R]; //vector of migrants
  vector[A*R] pop; //vector of log population
  int<lower=1> age[A*R]; //vector of ages 
  int<lower=1> region[A*R]; //vector of regions
}
parameters {
  vector[A] tau; //age-specific log rates
  vector[R] gamma; //region-specific intercept for log rates
  matrix[A,R] delta; //matrix of deviations
  real<lower=0> sigma; //sd of delta
}

transformed parameters{
  vector[N] log_rate; //vector of age & region specific log rate
  for (i in 1:N){
      log_rate[i] = tau[age[i]] + gamma[region[i]]+ delta[age[i],region[i]]; 
  }
  }


model {
  //likelihood
      y ~ poisson_log(log_rate + pop);

  //constraints
  for (a in 1:A){
    sum(delta[a, 1:R]) ~ normal(0, 0.01);
  }

  //priors
  tau ~ normal(0,1);
  gamma ~ normal(0,1);
  sigma ~ normal(0, 0.1);
  
  for (a in 1:A){
    for (r in 1:R){
      if (a==1){
        delta[a, r] ~ normal(0, sigma);
      } else{
        delta[a, r] ~ normal(delta[(a-1),r],sigma);
      }
    }
  }
}
