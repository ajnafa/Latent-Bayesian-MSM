// ARDL (t-1) Model with a Gaussian Likelihood for Simulation Comparison
data {
  int<lower = 0> N; // Observations
  vector[N] Y; // Outcome Stage Response
  int<lower = 1> K; // Number of Population-Level Effects and Intercept
  matrix[N, K] X; // Design Matrix for the Population-Level Effects
  
  // Coefficient Prior Scales
  vector[K] beta_sds;
}

transformed data {
  // Priors on the intercept and error term
  real alpha_prior_mu = mean(Y);
  real<lower = 0> alpha_prior_sd = 2 * sd(Y);
  real<lower = 0> sigma_prior = 1/sd(Y);
  
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
  vector[Kc] b;  // Population-Level Effects
  real Intercept;  // Intercept for the centered predictors
  real<lower=0> sigma;  // Dispersion parameter
}

model {
  // likelihood including constants
  target += normal_id_glm_lpdf(Y | Xc, Intercept, b, sigma);
  
  // Priors for the Model Parameters
  target += normal_lpdf(Intercept | alpha_prior_mu, alpha_prior_sd);
  target += normal_lpdf(b[1] | 0, beta_sds[2]);
  target += normal_lpdf(b[2] | 0, beta_sds[3]);
  target += normal_lpdf(b[3] | 0, beta_sds[4]);
  target += normal_lpdf(b[4] | 0, beta_sds[5]);
  target += normal_lpdf(b[5] | 0, beta_sds[6]);
  target += normal_lpdf(b[6] | 0, beta_sds[7]);
  target += exponential_lpdf(sigma | sigma_prior);
}

generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}

