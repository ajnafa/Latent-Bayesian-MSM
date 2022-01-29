/* Bayesian Single-Weighted Latent Inverse Probability Estimator with a 
* Gaussian Likelihood, Country Random Effects, and AR(1) errors
* Author: A. Jordan Nafa; Stan Version 2.28.2; Last Revised 1-28-2022 */
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
  
  // Calculate the mean for each observation-weight
  vector weights_mean(matrix weights, int N){
    vector[N] weights_mu;
    for (n in 1:N) {
      weights_mu[n] = mean(weights[, n]);
    }
    return weights_mu;
  }
  
  // Calculate the scale for each observation-weight
  vector weights_scale(matrix weights, int N){
    vector[N] weights_sd;
    for (n in 1:N) {
      weights_sd[n] = sd(weights[, n]);
    }
    return weights_sd;
  }
}

data {
  // Data for the Response Stage Model
  int<lower = 0> N; // Total Number of Observations
  vector[N] Y; // Response Vector
  int<lower = 1> K; // Number of Population-Level Effects
  matrix[N, K] X; // Design Matrix for the Population-Level Effects
  
  // Data for the Country level Effects
  int<lower = 1> J;           // Number of countries J
  int<lower = 1> jj[N];       // Observation-Country Mapping
  
  // Data from the Design Stage Model
  int<lower = 1> IPW_N; // Number of Rows in the Weights Matrix
  matrix[IPW_N, N] IPW; // Matrix of IP Weights from the Design Stage Model
  
  // Data for the ARMA Correlations
  int<lower=0> K_AR;  // AR order
  int<lower=0> J_lag[N]; // Number of Lags per Observation
}

transformed data {
  // Data for the Population-Level Latent Weights
  vector[N] ipw_mu; // Mean of the Population-Level IP Weights
  vector[N] ipw_sigma; // SD of the Population-Level IP Weights
  int max_lag = max(K_AR, 0);
  vector[N] Z_J = ones_vector(N); // Data for the Country level Intercepts
  
  // Calculate the location and scale for each observation weight
  ipw_mu = weights_mean(IPW, N);
  ipw_sigma = weights_scale(IPW, N);
}

parameters {
  vector[K] beta; // Population-Level Effects
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
  vector[N] mu = X * beta; 
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
  target += normal_lpdf(beta | 0, 2.5);
  target += normal_lpdf(ar | 0, 0.5) - 1 * log_diff_exp(normal_lcdf(1 | 0, 0.5), normal_lcdf(-1 | 0, 0.5));
  target += exponential_lpdf(sigma_upsilon | 0.5);
  target += normal_lpdf(z_upsilon | 0, 1);
  target += student_t_lpdf(sigma | 3, 0, 2.5) - 1 * student_t_lccdf(0 | 3, 0, 2.5);
}
