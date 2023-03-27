/* RW on overdispersion */
data {
  int<lower=0> A; //number of age groups
  int<lower=0> R; //number of regions
  int<lower=0> N; //number of data points
  int<lower=0> y[A*R]; //vector of migrants
  vector[A*R] pop; //vector of log population
  int<lower=1> age_vec[A*R];
  int<lower=1> region_vec[A*R];
  matrix[N,(A-1)] age; //matrix of ages 
  matrix[N,(R-1)] region; //matrix of regions
}
parameters {
  vector[A-1] tau; //age-specific log rates
  vector[R-1] gamma; //region-specific intercept for log rates
  real alpha;
  matrix[A,R] delta; //matrix of deviations
  real<lower=0> sigma; //sd of delta
}

transformed parameters{
  vector[N] log_rate; //vector of age & region specific log rate
  vector[N] fixed_eff;
  fixed_eff = alpha + age*tau + region*gamma; 
  for (i in 1:N){
      log_rate[i] = fixed_eff[i] + delta[age_vec[i],region_vec[i]] ; //age-specific rate + deviation
      }
    }
  
model {
  //likelihood
  y ~ poisson(exp(log_rate + pop));

  //constraints
  for (a in 1:A){
    sum(delta[a, 1:R]) ~ normal(0, 0.01);
  }

  //priors
  tau ~ normal(0,1);
  alpha ~ normal(0,1);
  gamma ~ normal(0,1);
  sigma ~ normal(0, 1);
  
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
