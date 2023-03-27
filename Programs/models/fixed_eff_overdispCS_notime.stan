/* Simple linear regression */
data {
  int<lower=0> A; //number of age groups
  int<lower=0> R; //number of regions
  int<lower=0> N; //number of data points
  int<lower=0> K; //number of knots
  int<lower=0> y[A*R]; //vector of migrants
  int<lower=0> pop[A*R]; //vector of population
  int<lower=1> age[A*R]; //vector of ages 
  int<lower=1> region[A*R]; //vector of regions
  matrix[A,K] B; 
}
parameters {
  vector[A] tau; //age-specific log rates
  vector[R] gamma; //region-specific intercept for log rates
  matrix[K,R] alpha; 
  real<lower=0> sigma_alpha; //sd of alpha
}

transformed parameters{
 //matrix[A,(R-1)] delta0;
 //matrix[A, 1] delta1;
 matrix[A,R] delta; //matrix of deviations
  matrix[A,R] log_rate; //vector of age & region specific log rate
  
  delta = B*alpha;
  //for (a in 1:A){
  //delta1[a,1] = -sum(delta0[a,1:(R-1)]);
  //}
  
  //delta = append_col(delta0, delta1) ;
  
  

  for (a in 1:A){
    for (r in 1:R){
      log_rate[a,r] = tau[a] + gamma[r] + delta[a,r]; //age-specific rate + deviation
    }
  }
}
model {
  //likelihood
    for (i in 1:N){
      y[i] ~ poisson_log(log_rate[age[i],region[i]] + log(pop[i]));
  }

  //constraints
  for (a in 1:A){
    sum(delta[a, 1:R]) ~ normal(0, 0.01);
  }
  
  for (r in 1:R){
    sum(delta[1:A, r]) ~ normal(0, 0.01);
  }

  //priors
  tau ~ normal(0,1);
  gamma ~ normal(0,1);
  sigma_alpha ~ normal(0,1);
  
  for (k in 1:K){
    for (r in 1:R){
      if (k==1){
        alpha[k, r] ~ normal(0, sigma_alpha);
      } else{
        alpha[k, r] ~ normal(alpha[(k-1),r],sigma_alpha);
      }
    }
  }
}