#----------Bayesian Latent Inverse Probability of Treatment Estimator----------#
#-Author: A. Jordan Nafa----------------------------Created: February 22, 2022-#
#-R Version: 4.1.2----------------------------------Revised: February 22, 2022-#

# Read in the weights matrix
ipw_matrix <- read_rds("data/iptw_weights_descs.rds") %>% 
  # Sort the data
  arrange(country, time)

# Weighted Log PDF of the Gaussian Pseudo-Likelihood----
gaussian_ipwt_stan_funs <- "
  // Weighted Log PDF of the Gaussian Pseudo-Likelihood
  real gaussian_ipw_lpdf(vector y, vector mu, real sigma, vector wts, int N) {
    real weighted_term;
    weighted_term = 0.00;
    for (n in 1:N) {
      weighted_term = weighted_term + wts[n] * (normal_lpdf(y[n] | mu[n], sigma));
    }
    return weighted_term;
  }
"

# Code for the latent inverse probability of treatment weights-----
ipwt_weights_mean <- "// Data from the Design Stage Model
vector[N] ipw_mu; // Mean of the Population-Level IP Weights"

ipwt_weights_scale <- 
  "vector[N] ipw_sigma; // Scale of the Population-Level IP Weights"

ipwt_weights_par <- 
  "vector<lower=0>[N] weights_z; // Standardized Latent IP Weights"

ipwt_weights_tpars <-
  "vector[N] w_tilde; // Latent Weights
  // Compute the Latent IPW Weights
  w_tilde = ipw_mu + ipw_sigma .* weights_z;"

# Add the Stan Code to the Model Block----
gaussian_ipwt_stan_vars <- stanvar(
  scode = gaussian_ipwt_stan_funs,
  block = "functions"
) + stanvar(
  x = ipw_matrix$mu_iptw,
  name = "ipw_mu",
  scode = ipwt_weights_mean,
  block = "data"
) + stanvar(
  x = ipw_matrix$sigma_iptw,
  name = "ipw_sigma",
  scode = ipwt_weights_scale,
  block = "data"
) + stanvar(
  scode = ipwt_weights_par,
  block = "parameters"
) + stanvar(
  scode = ipwt_weights_tpars,
  block = "tparameters"
)

# Define a function to calculate predictions----
posterior_predict_gaussian_ipw <- function(i, prep, ntrys = 5, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  sigma <- brms:::add_sigma_se(sigma, prep, i = i)
  ndraws <- prep$ndraws
  rcontinuous(
    n = ndraws, dist = "norm",
    mean = mu, sd = sigma,
    lb = prep$data$lb[i], ub = prep$data$ub[i],
    ntrys = ntrys
  )
}

# Define a function to calculate expectations----
posterior_epred_gaussian_ipw <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  mu
}

# Define a simple family for a gaussian marginal structural model
gaussian_ipw <- custom_family(
  name = "gaussian_ipw",
  dpars = c("mu", "sigma"),
  links = c("identity", "log"),
  type = "real",
  lb = c(NA, 0),
  vars = c("w_tilde", "N"),
  loop = FALSE
)

# Naive outcome Model, confounded response and unconfounded treatment
ipw_outcome_model_a <- bf(
  y_treated_conf ~ treat_bin_conf + beta + gamma + delta + 
    (1 | country) + ar(time = time, gr = country, p = 1),
  family = gaussian_ipw
)

# Priors for the model parameters (see prior_scaling function below)
outcome_priors_a <- prior(normal(52.09, 131.79441), class = "Intercept") +
  prior(normal(0, 4.11893), class = "b", coef = "treat_bin_conf") +
  prior(normal(0, 14.85291), class = "b", coef = "beta") +
  prior(normal(0, 8.34525), class = "b", coef = "gamma") +
  prior(normal(0, 14.70615), class = "b", coef = "delta") +
  prior(exponential(0.01518), class = "sd") +
  prior(normal(0, 0.5), class = "ar") +
  prior("target += exponential_lpdf(weights_z | 1);", check = FALSE)

# Test the model using brms----
ipwt_outcome_fit_a <- brm(
  ipw_outcome_model_a,
  data = sim_data_df,
  prior = outcome_priors_a,
  chains = 6, 
  cores = 6L,
  iter = 8000,
  seed = 12345, 
  stanvars = gaussian_ipwt_stan_vars,
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.9),
  refresh = 100,
  file = str_c(fits_dir, "IPWT_Gassian_AR"),
)
