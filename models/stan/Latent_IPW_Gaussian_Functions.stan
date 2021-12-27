/* Stan Functions for Bayesian Latent Inverse Probability Weighted Estimators 
* Author: A. Jordan Nafa; Stan Version 2.28.1; Last Revised 12-23-2021 */

// Weighted Log PDF of the Gaussian Pseudo-Likelihood
real normal_ipw_lpdf(vector y, vector mu, real sigma, vector w_tilde, int N) {
  real weighted_term;
  for (n in 1:N) {
    weighted_term = 0.00 + w_tilde[n] * normal_lpdf(y[n] | mu[n], sigma);
  }
  return weighted_term;
}
