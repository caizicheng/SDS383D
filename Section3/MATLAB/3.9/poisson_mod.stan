data {
  int N; // number of observations
  int P; // number of covariates
  matrix[N, P] X; //covariate matrix
  vector[N] y; //outcome vector
}

parameters {
  vector[P] beta; // the regression coefficients
}

model {
  // Define the priors
  beta ~ normal(0, 1);
  // The likelihood
  y ~ poisson(exp(X*beta));
}