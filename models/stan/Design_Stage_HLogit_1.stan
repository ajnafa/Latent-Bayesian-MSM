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
    X_sigma[k] = 1/sd(X[, k]);
  }
}

parameters {
  real alpha; // Population-Level Intercept
  vector[K] beta; // Regression Coefficients
  real<lower = 0> sigma; // Country-Level Standard Deviations
  vector[J] z; // Standardized Country-Level Effects
}

transformed parameters {
  vector[J] upsilon; // Actual Country-Level Effects
  upsilon = sigma * z;
}

model {
  // Priors for Model Parameters
  beta ~ normal(X_mu, X_sigma); // Priors for the coefficients
  sigma ~ exponential(0.5); // Prior for the country SDs
  z ~ std_normal(); // Prior for the standardized country effects
  alpha ~ student_t(4, 0, 1); // Prior on the population-level intercept
  
  // Bernoulli Likelihood
  target += bernoulli_logit_glm_lpmf(Y | X, alpha + upsilon[jj], beta);
}

generated quantities {
  vector[N] log_lik; // Pointwise Log Likelihood
  
  // Group Level Expectations
  real preds_groups[J] = bernoulli_logit_rng(alpha + upsilon[jj]);
  
  // Observation-Level Expectations
  real preds_obs[N] = bernoulli_logit_rng(alpha + upsilon[jj] + X * beta);

  // Calculate the Pointwise Log Likelihood
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(Y[n] | alpha + upsilon[jj[n]] + X[n] * beta);
  }
}
