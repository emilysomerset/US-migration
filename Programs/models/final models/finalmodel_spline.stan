/* RW on overdispersion */
data {
  int<lower=0> A; //number of age groups
  int<lower=0> R; //number of regions
  int<lower=0> N; //number of data points
  int<lower=0> T; //number of years 
  int<lower=0> K; //number of knots
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
  matrix[A,K] B; 
}
parameters {
  vector[A-1] tau; //age-specific log rates
  vector[R-1] gamma; //region-specific intercept for log rates
  real alpha;// intercept
  matrix[K,R] cs_coef; 
  matrix[T,R] delta2; //matrix of deviations
  vector[N] delta3; //matrix of deviations
  real<lower=0> sigma1; //sd of delta1
  real<lower=0> sigma2; //sd of delta2
}

transformed parameters{
  vector[N] log_rate; //vector of age & region specific log rate
  vector[N] fixed_eff;
  vector[N] deviation_fixed_eff;
  vector[N] both_deviation_fixed_eff;
  matrix[A,R] delta1; //matrix of deviations
  delta1 = B*cs_coef;
  fixed_eff = alpha + age*tau + region*gamma;
  for (i in 1:N){
      deviation_fixed_eff[i] = fixed_eff[i] + delta1[age_vec[i],region_vec[i]];
      both_deviation_fixed_eff[i] = fixed_eff[i] + delta1[age_vec[i],region_vec[i]]+ delta2[time_vec[i],region_vec[i]];
      log_rate[i] = fixed_eff[i] + delta1[age_vec[i],region_vec[i]] + delta2[time_vec[i],region_vec[i]] + delta3[i];
      }
    }
  

model {
  //likelihood
    target += poisson_lpmf(y| exp(log_rate+pop));
    target += normal_lpdf(tau|0,1);
    target += normal_lpdf(alpha|0,1);
    target += normal_lpdf(gamma|0,1);
    target += normal_lpdf(sigma1|0,1);
    for (k in 1:K){
      if (k==1){
        target += normal_lpdf(cs_coef[k,1:R]|0,sigma1);
      }else{
        target += normal_lpdf(cs_coef[k,1:R]|cs_coef[(k-1),1:R],sigma1);  //random walk on spline coefficients
    }}
    for (a in 1:A){
      target += normal_lpdf(sum(delta1[a,1:R])|0,0.01);  //constraint
    }
    target += normal_lpdf(sigma2|0,1);
    for (t in 1:T){
      if (t==1){
        target += normal_lpdf(delta2[t,1:R]|0,sigma2);
      }else{
        target += normal_lpdf(delta2[t, 1:R]|delta2[(t-1),1:R], sigma2);  //random walk on time
      }
    }
    // target += normal_lpdf(sigma3|0,1);
    target += normal_lpdf(delta3|0,0.3);
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
 pred_rate[i] = pred_fixed_eff[i]+ delta1[pred_age_vec[i],pred_region_vec[i]] + normal_rng(delta2[T,pred_region_vec[i]],sigma2) + normal_rng(0,0.3);
 }
}

