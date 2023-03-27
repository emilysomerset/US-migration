/* RW on overdispersion */
data {
  int<lower=0> A; //number of age groups
  int<lower=0> R; //number of regions
  int<lower=0> N; //number of data points
  int<lower=0> T; //number of years 
  int<lower=0> y[A*R*T]; //vector of migrants
  vector[A*R*T] pop; //vector of population
  int<lower=1> age[A*R*T]; //vector of ages 
  int<lower=1> region[A*R*T]; //vector of regions
  int<lower=1> time[A*R*T]; //vector of years
}
parameters {
  vector[A] tau; //age-specific log rates
  vector[R] alpha; //region-specific intercept for log rates
  vector[T] beta; //year-specific intercepts for log rates
  matrix[A,R] delta1; //matrix of deviations
  matrix[T,R] delta2; //matrix of deviations
  real<lower=0> sigma1; //sd of delta
  // real<lower=0> sigma2; //sd of delta
  real delta[A, T, R];
}

transformed parameters{
  vector[N] log_rate; //vector of age & region specific log rate
  for (i in 1:N){
      // log_rate[i] = tau[age[i]] + alpha[region[i]] + beta[time[i]] + delta1[age[i],region[i]] + delta2[time[i],region[i]]; //age-specific rate + deviation
      log_rate[i] = tau[age[i]] + alpha[region[i]] + beta[time[i]] + delta[age[i], time[i], region[i]]; //age-specific rate + deviation
      }
    }
  

model {
  //likelihood
    for (i in 1:N){
      y[i] ~ poisson_log(log_rate[i] + pop[i]);
  }

  //constraints
  // for (a in 1:A){
  //   sum(delta1[a, 1:R]) ~ normal(0, 0.01);
  // }
  //  for (t in 1:T){
  //   sum(delta2[t, 1:R]) ~ normal(0, 0.01);
  // }
  
  for (a in 1:A){
    for (t in 1:T){
      sum(delta[a,t,1:R]) ~ normal(0, 0.01);
    }
  }

  //priors
  tau ~ normal(0,1);
  alpha ~ normal(0,1);
  beta ~ normal(0,1);
  sigma1 ~ normal(0, 1);
  // sigma2 ~ normal(0, 1);

  
  // for (a in 1:A){
  //   for (r in 1:R){
  //     if (a==1){
  //       delta1[a, r] ~ normal(0, sigma1);
  //     } else{
  //       delta1[a, r] ~ normal(delta1[(a-1),r],sigma1);
  //     }
  //   }
  // }
  // 
  //   for (t in 1:T){
  //   for (r in 1:R){
  //     if (t==1){
  //       delta2[t, r] ~ normal(0, sigma2);
  //     } else{
  //       delta2[t, r] ~ normal(delta2[(t-1),r],sigma2);
  //     }
  //   }
  // }
  
  for (a in 1:A){
    for (t in 1:T){
      for (r in 1:R){
      if (a==1 & t==1){
        delta[a,t,r] ~ normal(0, sigma1);
      } else{
        delta[a,t,r] ~ normal(delta[(a-1),(t-1),r])
      }
      }
    }
  }
}
