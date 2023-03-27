/* RW on overdispersion */
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
  matrix[N,(A-1)*(R-1)] age_region; //interaction
  int<lower=1> age_vec[A*R*T];
  int<lower=1> time_vec[A*R*T]; 
  int<lower=1> region_vec[A*R*T]; 
  matrix[A*R,(A-1)] pred_age; //matrix of ages 
  matrix[A*R,(R-1)] pred_region; //matrix of regions
  matrix[A*R,(A-1)*(R-1)] pred_age_region; //interaction
  int<lower=1> pred_age_vec[A*R];
  int<lower=1> pred_region_vec[A*R];
}
parameters {
  vector[A-1] tau; //age-specific log rates
  vector[R-1] gamma; //region-specific intercept for log rates
  vector[(A-1)*(R-1)] xi;
  real alpha;// intercept
  matrix[T,A] delta1; //matrix of deviations
  matrix[T,R] delta2; //matrix of deviations
  // matrix[A,R] delta3; //matrix of deviations
  real<lower=0> sigma1; //sd of delta1
  real<lower=0> sigma2; //sd of delta2
  // real<lower=0> sigma3;
}

transformed parameters{
  vector[N] log_rate; //vector of age & region specific log rate
  vector[N] fixed_eff;
  fixed_eff = alpha + age*tau + region*gamma + age_region*xi;
  for (i in 1:N){
      log_rate[i] = fixed_eff[i] + delta1[time_vec[i],age_vec[i]] + delta2[time_vec[i],region_vec[i]];
      }
    }
  

model {
  //likelihood
      y ~ poisson_log(log_rate + pop);

  // //constraints
  // for (t in 1:T){
  //   sum(delta1[t, 1:A]) ~ normal(0, 0.001);
  //   sum(delta2[t, 1:R]) ~ normal(0, 0.001);
  // }
  
  // for (a in 1:A){
  //   sum(delta3[a, 1:R]) ~ normal(0, 0.001);
  // }


  //priors
  tau ~ normal(0,1);
  alpha ~ normal(0,1);
  xi ~ normal(0,1);
  sigma1 ~ normal(0, 1);
  sigma2 ~ normal(0, 1);

  
delta1[1,1:A] ~ normal(0, sigma1);
delta2[1,1:R] ~ normal(0, sigma2);
sum(delta1[1, 1:A]) ~ normal(0, 0.001);
sum(delta2[1, 1:R]) ~ normal(0, 0.001);

for (t in 2:T){
    sum(delta1[t, 1:A]) ~ normal(0, 0.001);
    sum(delta2[t, 1:R]) ~ normal(0, 0.001);
    delta1[t, 1:A] ~ normal(delta1[(t-1), 1:A], sigma1);
    delta2[t, 1:R] ~ normal(delta2[(t-1), 1:R], sigma2);
}

  
  // for (a in 1:A){
  //   for (r in 1:R){
  //     if (a==1){
  //       delta3[a, r] ~ normal(0, sigma3);
  //     } else{
  //       delta3[a, r] ~ normal(delta3[(a-1),r],sigma3);
  //     }
  //   }
  // }
}
  
  
generated quantities {

vector[N] rate_p;
vector[A*R] pred_rate;
vector[A*R] pred_fixed_eff;

for (i in 1:N){
if (time_vec[i]==1){rate_p[i] = 0;
} else{
  rate_p[i] = fixed_eff[i] + normal_rng(delta1[(time_vec[i]-1),age_vec[i]],sigma1)+ normal_rng(delta2[(time_vec[i]-1),region_vec[i]],sigma2);
}}

pred_fixed_eff = alpha + pred_age*tau + pred_region*gamma+ pred_age_region*xi;
for (i in 1:(A*R)){
 pred_rate[i] = pred_fixed_eff[i]+ normal_rng(delta1[T,pred_age_vec[i]],sigma1) + normal_rng(delta2[T,pred_region_vec[i]],sigma2);
 }
}

