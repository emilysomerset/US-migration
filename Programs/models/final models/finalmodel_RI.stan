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
  int<lower=1> age_vec[A*R*T];
  int<lower=1> time_vec[A*R*T]; 
  int<lower=1> region_vec[A*R*T]; 
  matrix[A*R,(A-1)] pred_age; //matrix of ages 
  matrix[A*R,(R-1)] pred_region; //matrix of regions
  int<lower=1> pred_age_vec[A*R];
  int<lower=1> pred_region_vec[A*R];
}
parameters {
  vector[A-1] tau; //age-specific log rates
  vector[R-1] gamma; //region-specific intercept for log rates
  real alpha;// intercept
  vector[N] delta3; //matrix of deviations
  real<lower=0> sigma3; //sd of delta2
}

transformed parameters{
  vector[N] log_rate; //vector of age & region specific log rate
  vector[N] fixed_eff;
  fixed_eff = alpha + age*tau + region*gamma;
  for (i in 1:N){
      log_rate[i] = fixed_eff[i] + delta3[i];
      }
    }
  

model {
  //likelihood
    target += poisson_lpmf(y| exp(log_rate+pop));
    target += normal_lpdf(tau|0,1);
    target += normal_lpdf(alpha|0,1);
    target += normal_lpdf(gamma|0,1);
    target += normal_lpdf(sigma3|0,1);
    target += normal_lpdf(delta3|0,sigma3);
}
  
  
generated quantities {

// vector[N] rate_p;
int<lower=0> y_p[N];
vector[N] generated_rate;
vector[A*R] pred_rate;
vector[A*R] pred_fixed_eff;
real sum_delta3;
sum_delta3 = sum(delta3);

y_p = poisson_rng(exp(log_rate + pop));

for (i in 1:N){
generated_rate[i] = y_p[i]/exp(pop[i]);
}

// predicting year 2020
pred_fixed_eff = alpha + pred_age*tau + pred_region*gamma;
for (i in 1:(A*R)){
 pred_rate[i] = pred_fixed_eff[i]+ delta1[pred_age_vec[i],pred_region_vec[i]] + normal_rng(delta2[T,pred_region_vec[i]],sigma2[pred_region_vec[i]]) + normal_rng(0,sigma3);
 }
}

