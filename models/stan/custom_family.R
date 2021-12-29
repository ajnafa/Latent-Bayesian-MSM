#-------------Custom brms Families for a Pseudo-Bayesian IPW Estimator----------
#-Author: A. Jordan Nafa----------------------------Created: December 28, 2021-#
#-R Version: 4.1.0----------------------------------Revised: December 28, 2021-#

# Load the necessary libraries
library(brms) # Version 2.16.4

#------------------------------------------------------------------------------#
#--------------------------Stan Functions and Code------------------------------
#------------------------------------------------------------------------------#

# Stan Code for the weighted likelihood
"// Weighted Log PDF of the Gaussian Pseudo-Likelihood
real normal_ipw_lpdf(vector y, vector mu, real sigma, vector w_tilde, int N) {
  real weighted_term;
  weighted_term = 0.00;
  for (n in 1:N) {
    weighted_term = weighted_term + w_tilde[n] * (normal_lpdf(y[n] | mu[n], sigma));
  }
  return weighted_term;
}" -> stan_funs

# Add the Stan Code to the functions block
stan_funs_vars <- stanvar(
  scode = stan_funs, 
  block = "functions"
)

# Add the Stan Code to the functions block
stan_data_vars_ipwn <- stanvar(
  x = 8000, # Set this to the number of rows in the weights matrix
  name = "IPW_N",
  scode = "int<lower = 1> IPW_N; // Number of Rows in the Weights Matrix", 
  block = "data"
)

stan_data_vars_ipwmat <- stanvar(
  x = ipw_mat, # An 8000 x N matrix of IP Weights
  name = "IPW",
  scode = "matrix[IPW_N, N] IPW; // Matrix of IP Weights from the Design Stage Model", 
  block = "data"
)

# Defining the IPW Matrix in the transformed data block
tdata_vars <- stanvar(
  scode = "// Data for the Population-Level Latent Weights
  vector[N] weights_mean; // Mean of the Population-Level IP Weights
  vector[N] weights_sd; // SD of the Population-Level IP Weights
  // Calculate the location and scale for each observation weight
  for (n in 1:N) {
    weights_mean[n] = mean(IPW[, n]);
    weights_sd[n] = sd(IPW[, n]);
  }",
  block = "tdata"
)

# Defining the latent weights as a transformed parameter
pars_vars <- stanvar(
  scode = "vector<lower=0>[N] weights_z; // Standardized Latent IP Weights",
  block = "parameters"
)

tpars_vars <- stanvar(
  scode = "vector[N] w_tilde; // Latent IPW Weights
  // Compute the Latent IPW Weights
  w_tilde = weights_mean + weights_sd .* weights_z;",
  block = "tparameters"
)

# Combine all of these in a single object
normal_ipw_stanfuns <- c(
  stan_funs_vars,
  stan_data_vars_ipwn,
  stan_data_vars_ipwmat,
  tdata_vars,
  pars_vars,
  tpars_vars
)


# Define a simple family for a marginal structural model
normal_ipw <- custom_family(
  name = "normal_ipw",
  dpars = c("mu", "sigma"),
  links = c("identity", "log"),
  type = "real",
  lb = c(NA, 0),
  vars = c("w_tilde", "N"),
  loop = FALSE
)

# @todo: Generalize to multilevel setting, write up additional functions
# for predictions, expectations, ect; test on simulated data with varying
# effect size

#------------------------------------------------------------------------------#
#----------------------------Custom Family Example------------------------------
#------------------------------------------------------------------------------#


# Priors for the model
normal_ipw_priors <- 
  prior("target += exponential_lpdf(weights_z | 1)", check = FALSE) +
  prior("normal(0, 5)", class = "b") +
  prior("student_t(3, 0, 14.8)", class = "sigma") +
  prior("student_t(3, 31, 14.8)", class = "Intercept")


# Fit the model using brms----
outcome_model <- brm(
  malaria_risk ~ net_num,
  family = normal_ipw,
  data = nets,
  prior = normal_ipw_priors,
  chains = 6, 
  cores = 6L,
  data2 = list(IPW_N = 8000, IPW = ipw_mat),
  iter = 8000, # Takes about 140 seconds on my windows laptop
  seed = 1234, 
  backend = "cmdstanr",
  stanvars = normal_ipw_stanfuns,
  save_pars = save_pars(all = TRUE),
  refresh = 100
)
