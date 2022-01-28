/* Design/Treatment Stage Multilevel Logitistic Regression Model 1
Author: A. Jordan Nafa; Last Revised 01/27/2020 under Stan 2.8.2*/
data {
  int<lower = 1> N;           // Observations N
  int Y[N];                 // Response Vector
  int<lower = 1> K;           // Number of Fixed Effects
  matrix[N, K] X;           // Centered Design Matrix for the Fixed Effects
  // Data for the Country level Effects
  int<lower = 1> J;           // Number of countries J
  int<lower = 1> jj[N];       // Observation-Country Mapping
}

transformed data {
  vector[K] X_sigma;          // SDs for Coefficient Priors
  vector[K] X_mu = zeros_vector(K); // Means for Coefficient Priors
  vector[N] Z_J = ones_vector(N); // Data for the Country level Intercepts
  
  // Calculate the standard deviation of each column in X
  for (k in 1:K) {
    X_sigma[k] = sd(X[, k]);
  }
}

parameters {
  real alpha; // Population-Level Intercept
  vector[K] beta; // Regression Coefficients
  real<lower = 0> sigma_upsilon; // Country-Level Standard Deviations
  vector[J] z_upsilon; // Standardized Country-Level Effects
}

transformed parameters {
  vector[J] upsilon = sigma_upsilon * z_upsilon; // Actual Country-Level Effects
  vector[N] mu; // Linear Predictor
  mu = alpha + upsilon[jj] .* Z_J;
}

model {
  // Priors for Model Parameters
  target += exponential_lpdf(sigma_upsilon | 1); // Prior for the country SDs
  target += normal_lpdf(z_upsilon | 0, 1); // Prior for the standardized country effects
  target += normal_lpdf(alpha | 0, 3); // Prior on the population-level intercept
  target += normal_lpdf(beta | X_mu, X_sigma); // Priors for the coefficients
  
  // Bernoulli Likelihood
  target += bernoulli_logit_glm_lpmf(Y | X, mu, beta);
}

generated quantities {
  vector[N] log_lik; // Pointwise Log Likelihood
  vector[N] treat_preds; // Observation-Level Expectations
  
  for(n in 1:N) {
    treat_preds[n] = bernoulli_logit_rng(mu[n] + X[n] * beta);
  }
  
  // Calculate the Pointwise Log Likelihood
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(Y[n] | alpha + X[n] * beta);
  }
}
