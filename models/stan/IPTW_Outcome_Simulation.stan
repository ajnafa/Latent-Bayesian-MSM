/* Pseudo Bayesian Inverse Probability of Treatment Weighted Estimator
* Author: A. Jordan Nafa; Stan Version 2.30.1; Last Revised 08-29-2022 */
functions {
  // Weighted Log PDF of the Gaussian Pseudo-Likelihood
  real normal_ipw_lpdf(vector y, vector mu, real sigma, vector w_tilde, int N) {
    real weighted_term;
    weighted_term = 0.00;
    for (n in 1:N) {
      weighted_term = weighted_term + w_tilde[n] * (normal_lpdf(y[n] | mu[n], sigma));
    }
    return weighted_term;
  }
}

data {
  int<lower = 0> N; // Observations
  vector[N] Y; // Outcome Stage Response
  int<lower = 1> K; // Number of Population-Level Effects and Intercept
  matrix[N, K] X; // Design Matrix for the Population-Level Effects
  
  // Statistics from the Design Stage Model
  vector<lower = 0>[N] ipw_mu; // Mean of the Population-Level Weights
  vector<lower = 0>[N] ipw_sigma; // Scale of the Population-Level Weights
  
  // Prior on the scale of the weights
  real<lower = 0> sd_prior_shape1;
  real<lower = 0> sd_prior_shape2;
  
  // Prior on the coefficients and intercept
  real<lower = 0> b_prior_sigma;
  real alpha_prior_mu;
  real<lower = 0> alpha_prior_sigma;
  real<lower = 0> sigma_prior;
}

transformed data {
  int Kc = K - 1;   
  matrix[N, Kc] Xc;  // Centered version of X without an Intercept
  vector[Kc] means_X;  // Column Means of the Uncentered Design Matrix
  
  // Centering the design matrix
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}

parameters {
  vector[Kc] b; // Population-Level Effects
  real Intercept; // Population-Level Intercept for the Centered Predictors
  real<lower=0> sigma;  // Dispersion Parameter
  real<lower=0, upper=1> weights_z[N]; // Parameter for the IPT Weights
}

transformed parameters {
  // Compute the IPT Weights
  vector[N] w_tilde; // IPT Weights
  w_tilde = ipw_mu + ipw_sigma * weights_z[1];
}

model {
  // Initialize the Linear Predictor
  vector[N] mu = Intercept + Xc * b;
  
  // Sampling the Weights
  target += beta_lpdf(weights_z | sd_prior_shape1, sd_prior_shape2);
  
  // Priors for the Model Parameters
  target += normal_lpdf(Intercept | alpha_prior_mu, alpha_prior_sigma);
  target += normal_lpdf(b | 0, b_prior_sigma);
  target += exponential_lpdf(sigma | sigma_prior);
  
  // Weighted Likelihood
  target += normal_ipw_lpdf(Y | mu, sigma, w_tilde, N);
}

generated quantities {
  // Population-level Intercept on the Original Scale
  real b_Intercept = Intercept - dot_product(means_X, b);
}

  
