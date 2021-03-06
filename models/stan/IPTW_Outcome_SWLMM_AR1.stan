/* Bayesian Single-Weighted Latent Inverse Probability Estimator with a 
* Gaussian Likelihood, Country Random Effects, and AR(1) errors
* Author: A. Jordan Nafa; Stan Version 2.28.2; Last Revised 1-29-2022 */
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
  // Data for the Response Stage Model
  int<lower = 0> N; // Total Number of Observations
  vector[N] Y; // Response Vector
  int<lower = 1> K; // Number of Population-Level Effects
  matrix[N, K] X; // Design Matrix for the Population-Level Effects
  
  // Data for the Country level Effects
  int<lower = 1> J;   // Number of Countries J
  int<lower = 1> jj[N]; // Observation-Country Mapping
  
  // Data from the Design Stage Model
  vector[N] ipw_mu; // Mean of the Population-Level IP Weights
  vector[N] ipw_sigma; // Scale of the Population-Level IP Weights
  
  // Data for the ARMA Correlations
  int<lower = 0> K_AR;  // AR Order
  int<lower = 0> K_MA; // MA Order
  int<lower = 0> J_lag[N]; // Number of Lags per Observation
}

transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // Centered version of X without an Intercept
  vector[Kc] means_X;  // Column Means of the Uncentered Design Matrix
  
  int max_lag = max(K_AR, K_MA);
  vector[N] Z_J = ones_vector(N); // Data for the Country level Intercepts
  
  // Centering the Design Matrix
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}

parameters {
  vector[Kc] beta; // Population-Level Effects
  real alpha; // Population-Level Intercept for the Centered Predictors
  real<lower=0> sigma;  // Dispersion Parameter
  vector<lower=0>[N] weights_z; // Standardized Latent IP Weights
  vector[K_AR] ar;  // Autoregressive Coefficients
  real<lower = 0> sigma_upsilon; // Country-Level Standard Deviations
  vector[J] z_upsilon; // Standardized Country-Level Effects
}

transformed parameters {
  vector[J] upsilon = sigma_upsilon * z_upsilon; // Actual Country-Level Effects
  vector[N] w_tilde; // Latent IPT Weights
  
  // Compute the Latent IPW Weights
  w_tilde = ipw_mu + ipw_sigma .* weights_z;
}

model {
  // Initialize a Matrix for Storing the Lagged Residuals
  matrix[N, max_lag] Err = rep_matrix(0, N, max_lag);
  vector[N] err;  // Actual Residuals
  
  // Initialize the Linear Predictor
  vector[N] mu = alpha + Xc * beta; 
  profile ("Linear Predictor") {
    for (n in 1:N) {
      // Add Terms to the linear predictor
      mu[n] += upsilon[jj[n]] * Z_J[n];
    }
  }
  // Include AR(1) Terms
  profile ("AR(1) Terms") {
    for (n in 1:N) {
      err[n] = Y[n] - mu[n];
      for (i in 1:J_lag[n]) {
        Err[n + 1, i] = err[n + 1 - i];
      }
      mu[n] += Err[n, 1:K_AR] * ar;
    }
  }
  // Likelihood
  profile ("Weighted Likelihood") {
    target += normal_ipw_lpdf(Y | mu, sigma, w_tilde, N);
  }
  // Sampling the Weights
  target += exponential_lpdf(weights_z | 1);
  
  // Priors for the Model Parameters
  target += normal_lpdf(alpha | 54, 131);
  target += normal_lpdf(beta[1] | 0, 4);
  target += normal_lpdf(beta[2] | 0, 15);
  target += normal_lpdf(beta[3] | 0, 8);
  target += normal_lpdf(beta[4] | 0, 15);
  target += normal_lpdf(ar | 0, 0.5) 
  - 1 * log_diff_exp(normal_lcdf(1 | 0, 0.5), normal_lcdf(-1 | 0, 0.5));
  target += exponential_lpdf(sigma_upsilon | 0.015);
  target += normal_lpdf(z_upsilon | 0, 1);
  target += exponential_lpdf(sigma | 0.015);
}

generated quantities {
  // Population-level Intercept on the Original Scale
  real Intercept = alpha - dot_product(means_X, beta);
}
