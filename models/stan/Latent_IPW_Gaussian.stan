/* Bayesian Latent Inverse Probability Estimator with a Gaussian Likelihood
* Author: A. Jordan Nafa; Stan Version 2.28.1; Last Revised 12-23-2021 */
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
  // Data for the outcome model
  int<lower = 0> N; // Total Number of Observations
  vector[N] Y; // Response Vector
  int<lower = 1> K; // Number of Population-Level Effects
  matrix[N, K] X; // Design Matrix for the Population-Level Effects
  // Data from the Design Stage Model
  int<lower = 1> IPW_N; // Number of Rows in the Weights Matrix
  matrix[IPW_N, N] IPW; // Matrix of IP Weights from the Design Stage Model
}
transformed data {
  // Data for the Population-Level Latent Weights
  vector[N] ipw_mu; // Mean of the Population-Level IP Weights
  vector[N] ipw_sigma; // SD of the Population-Level IP Weights
  // Calculate the location and scale for each observation weight
  for (n in 1:N) {
    ipw_mu[n] = mean(IPW[, n]);
    ipw_sigma[n] = sd(IPW[, n]);
  }
  // Centering the Predictor Matrix
  int Kc = K - 1;
  matrix[N, Kc] Xc; // Centered version of X without an Intercept
  vector[Kc] means_X; // Column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  real alpha; // Temporary Population-Level Intercept
  vector[Kc] beta; // Population-Level Effects
  real<lower=0> sigma;  // Dispersion Parameter
  vector<lower=0>[N] weights_z; // Standardized Latent IP Weights
}
transformed parameters {
  vector[N] w_tilde; // Latent IPW Weights
  // Compute the Latent IPW Weights
  w_tilde = ipw_mu + ipw_sigma .* weights_z;
}
model {
  // Likelihood
  vector[N] mu = alpha + Xc * beta;
  target += normal_ipw_lpdf(Y | mu, sigma, w_tilde, N);
  // Sampling the Weights
  target += exponential_lpdf(weights_z | 1);
  // Priors for the model parameters
  target += student_t_lpdf(alpha | 3, 0, 2.5);
  target += normal_lpdf(beta | 0, 5);
  target += student_t_lpdf(sigma | 3, 0, 10) - 1 * student_t_lccdf(0 | 3, 0, 10);
}
generated quantities {
  // Actual Population-Level Intercept
  real Intercept = alpha - dot_product(means_X, beta);
  // Calculating Treatment Effects
  real yhat_treated = Intercept + beta[1]*1; // Predictions for the Treated Units
  real yhat_untreated = Intercept + beta[1]*0; // Predictions for the Untreated Units
  real mu_ate = yhat_treated - yhat_untreated; // Average difference between groups
  // Pseudo Log-Likelihood
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = w_tilde[n] * (normal_lpdf(Y[n] | alpha + Xc[n] * beta, sigma));
  }
}
