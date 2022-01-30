#-------------------Gaussian Marginal Structural Models-------------------------
#-Author: A. Jordan Nafa-----------------------------Created: January 28, 2022-#
#-R Version: 4.1.0-----------------------------------Revised: January 28, 2022-#

# Set Project Options----
options(
  digits = 4, # Significant figures output
  scipen = 999, # Disable scientific notation
  repos = getOption("repos")["CRAN"],
  mc.cores = 12
)

# Load the necessary libraries----
pacman::p_load(
  "tidyverse",
  "data.table",
  "dtplyr",
  "brms",
  "tidybayes",
  "cmdstanr"
)

# Load the helper functions
source("scripts/00_Helper_Functions.R")

# Specify the directory for stan models
stan_dir <- "models/stan/"

# Specify the directory for fitted model files
fits_dir <- "models/fits/outcome_stage/"

#------------------------------------------------------------------------------#
#------------------------------Data Preparation---------------------------------
#------------------------------------------------------------------------------#

# Load the simulated data
sim_data_ls <- read_rds("data/gaussian_ar1_sims.rds")

# Take a look at the loaded data structure
glimpse(sim_data_ls$data)

# Get the true values for the treatment
(true_treat <- sim_data_ls$true_treat)
# intercept   beta     gamma     delta 
# 0.59736   1.80015  -0.03517  -0.04994 

# Get the true values for the treatment
(true_resp <- sim_data_ls$true_resp)
# intercept     treat      beta     gamma     delta 
#  12.359       6.856     3.188    -2.449     5.406 

# Check distribution of the treatment
xtabs(~ country + treat_bin_conf, sim_data_ls$data)

# Extract the data
sim_data_df <- sim_data_ls$data %>%
  # Group the data by country
  group_by(country) %>%
  # Centering and scaling the predictors
  mutate(
    across(
      beta:delta,
      list(
        wi = ~ .x - mean(.x),
        be = ~ mean(.x)
      )
    )) %>% 
  # Ungroup the data
  ungroup() %>% 
  # Centering the scaling group-level predictors
  mutate(country = as_factor(country)) %>% 
  # Sort the data
  arrange(country, time)
  

# Read in the weights matrix
ipw_descs <- read_rds("data/iptw_weights_descs.rds") %>% 
  # Sort the data
  arrange(country, time)

#------------------------------------------------------------------------------#
#-------------Estimate the Naive Outcome Stage Model for URTC-------------------
#------------------------------------------------------------------------------#

# Naive outcome Model, confounded response and unconfounded treatment
naive_outcome_model_a <- bf(
  y_treated_conf ~ treat_bin_conf + beta + gamma + delta + 
    (1 | country) + ar(time = time, gr = country, p = 1),
  family = gaussian(link = "identity")
)

# Priors for the model parameters (see prior_scaling function below)
outcome_priors_a <- prior(normal(52.09, 131.79441), class = "Intercept") +
  prior(normal(0, 4.11893), class = "b", coef = "treat_bin_conf") +
  prior(normal(0, 14.85291), class = "b", coef = "beta") +
  prior(normal(0, 8.34525), class = "b", coef = "gamma") +
  prior(normal(0, 14.70615), class = "b", coef = "delta") +
  prior(exponential(0.01518), class = "sd") +
  prior(normal(0, 0.5), class = "ar")

# Fit the model using brms----
naive_outcome_fit_a <- brm(
  naive_outcome_model_a,
  data = sim_data_df,
  prior = outcome_priors_a,
  chains = 6, 
  cores = 6L,
  iter = 8000,
  seed = 12345, 
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.9),
  refresh = 100,
  save_model = str_c(stan_dir, "Naive_Outcome_LMM.stan"),
  file = str_c(fits_dir, "URTC/Naive_Outcome_LMM_URTC"),
)

#------------------------------------------------------------------------------#
#-------------Estimate the Naive Outcome Stage Model for CRUT-------------------
#------------------------------------------------------------------------------#

# Naive Outcome Model, confounded response and unconfounded treatment
naive_outcome_model_b <- bf(
  y_conf_treated ~ treat_bin + beta + gamma + delta + 
    (1 | country) + ar(time = time, gr = country, p = 1),
  family = gaussian(link = "identity")
)

# Priors for the model parameters (see prior_scaling function below)
outcome_priors_b <- prior(normal(54.17, 132.4), class = "Intercept") +
  prior(normal(0, 4.11893), class = "b", coef = "treat_bin") +
  prior(normal(0, 14.85291), class = "b", coef = "beta") +
  prior(normal(0, 8.34525), class = "b", coef = "gamma") +
  prior(normal(0, 14.70615), class = "b", coef = "delta") +
  prior(exponential(0.0151), class = "sd") +
  prior(normal(0, 0.5), class = "ar")

# Fit the model using brms----
naive_outcome_fit_b <- brm(
  naive_outcome_model_b,
  data = sim_data_df,
  prior = outcome_priors_b,
  chains = 6, 
  cores = 6L,
  iter = 8000,
  seed = 12345, 
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.9),
  refresh = 100,
  save_model = str_c(stan_dir, "Naive_Outcome_LMM_CRUT.stan"),
  file = str_c(fits_dir, "CRUT/Naive_Outcome_LMM_CRUT")
)

#------------------------------------------------------------------------------#
#-------------Estimate the Naive Outcome Stage Model for CRCT-------------------
#------------------------------------------------------------------------------#

# Naive Outcome Model, confounded response and confounded treatment
naive_outcome_model_c <- bf(
  y_conf_treated_conf ~ treat_bin_conf + beta + gamma + delta + (1 | country),
  family = gaussian(link = "identity")
)

# Priors for the model parameters (see prior_scaling function below)
outcome_priors_c <- prior(normal(0, 132.4), class = "Intercept") +
  prior(normal(0, 4), class = "b", coef = "treat_bin_conf") +
  prior(normal(0, 14), class = "b", coef = "beta") +
  prior(normal(0, 8), class = "b", coef = "gamma") +
  prior(normal(0, 14), class = "b", coef = "delta") +
  prior(normal(0, 132.4), class = "sd")

# Fit the model using brms----
naive_outcome_fit_c <- brm(
  naive_outcome_model_c,
  data = sim_data_df,
  prior = outcome_priors_c,
  chains = 6, 
  cores = 6L,
  iter = 8000,
  seed = 12345, 
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  control = list(adapt_delta = 0.9),
  refresh = 100,
  save_model = str_c(stan_dir, "Naive_Outcome_LMM_CRCT.stan"),
  file = str_c(fits_dir, "CRCT/Naive_Outcome_LMM_CRCT_3")
)

#------------------------------------------------------------------------------#
#----------Estimate the Latent IPTW Single-Weighted Outcome Stage Model---------
#------------------------------------------------------------------------------#

# Take advantage of brms functionality because I'm lazy
mod_data <- make_standata(
  y_conf_treated_conf ~ treat_bin_conf + beta + gamma + delta +
    (1 | country) + ar(time = time, gr = country, p = 1), 
  data = sim_data_df
  )

# Get priors for the model parameters
prior_scaling(y = mod_data$Y, df = mod_data$X)

# Prepare the data for use with Stan
iptw_stan_data <- list(
  N = nrow(sim_data_df),
  K = mod_data$K,
  Y = mod_data$Y,
  X = mod_data$X,
  J = mod_data$N_1,
  jj = mod_data$J_1,
  K_AR = 1,
  K_MA = 0,
  J_lag = mod_data$J_lag,
  ipw_mu = ipw_descs$mu_iptw,
  ipw_sigma = ipw_descs$sigma_iptw
)

# Load the Stan Model File for the Single Weighted Model
SWLMM_AR1 <- str_c(stan_dir, "IPTW_Outcome_SWLMM_AR1.stan")

# Compile the Stan model
swlmm_ar1_mod <- cmdstan_model(
  SWLMM_AR1, 
  dir = str_c(fits_dir, "SWLMM AR1/"),
  force_recompile = TRUE
)

# Print the model code
str_split(swlmm_ar1_mod$code(), pattern = ";", simplify = T)

# Fit the Outcome-Stage Model; Run time is approximately 10 minutes
swlmm_ar1_fit <- swlmm_ar1_mod$sample(
  data = iptw_stan_data,
  seed = 123456,
  refresh = 50,
  output_dir = str_c(fits_dir, "SWLMM AR1/"),
  sig_figs = 5,
  parallel_chains = 6,
  chains = 6,
  iter_warmup = 5000,
  iter_sampling = 5000,
  max_treedepth = 11
)

# Write the model object to an RDS file
swlmm_ar1_fit$save_object(
  file = str_c(fits_dir, "SWLMM AR1/IPTW_Outcome_SWLMM_AR1_CRCT.rds")
  )

# Extract the Posterior Draws
draws <- swlmm_ar1_fit$draws(
  variables = c("Intercept", "beta", "upsilon", "ar"),
  format = "df"
  )

# Print a Summary of the Posterior Draws
(summ_draws <- draws %>% summarise_draws(.cores = 12))
