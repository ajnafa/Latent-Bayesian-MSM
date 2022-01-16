/* Stan Functions for Bayesian Latent Inverse Probability Weighted Estimators 
* Author: A. Jordan Nafa; Stan Version 2.28.1; Last Revised 12-23-2021 */

// Weighted Log PDF of the Gaussian Pseudo-Likelihood
real normal_ipw_lpdf(vector y, vector mu, real sigma, vector w_tilde, int N) {
  real weighted_term;
  weighted_term = 0.00;
  for (n in 1:N) {
    weighted_term = weighted_term + w_tilde_i[n] * (normal_lpdf(y[n] | mu[n], sigma));
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
