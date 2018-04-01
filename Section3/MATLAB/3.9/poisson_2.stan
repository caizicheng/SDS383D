 // Basic Poisson glm
 
 data {
   // Define variables in data
   // Number of observations (an integer)
   int<lower=0> N;

   
 
   // Covariates
   int <lower=0, upper=1> intercept[N];
   int <lower=-1, upper=12> grade[N];
   int <lower=0, upper=1> se_attend[N];
   int <lower=0, upper=1> gender[N];
   int <lower=0, upper=1> genderattendance[N];
   int <lower=-1, upper=12> gendergrade[N];
   // Count outcome
   int<lower=5> y[N];

 }
 
 parameters {
   // Define parameters to estimate
   real beta[6];
 }
 
 transformed parameters  {
   //
   real lp[N];
   real <lower=0> mu[N];
 
   for (i in 1:N) {
     // Linear predictor
     lp[i] = beta[1] + beta[2]*grade[i] + beta[3]*gender[i] + beta[4]*se_attend[i]+beta[5]*genderattendance[i]+beta[6]*gendergrade[i];
 
     // Mean
     mu[i] = exp(lp[i]);
   }
 }
 
 model {
   // Prior part of Bayesian inference
   beta[1]~normal(0,1);
   beta[2]~normal(0,1);
   beta[3]~normal(0,1);
   beta[4]~normal(0,1);
   beta[5]~normal(0,1);
   beta[6]~normal(0,1);
 
   // Likelihood part of Bayesian inference
   y ~ poisson(mu);
 }