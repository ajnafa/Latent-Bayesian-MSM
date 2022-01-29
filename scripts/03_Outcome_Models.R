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
      c(beta:delta),
      list(
        wi = ~ (.x - mean(.x))/(2*sd(.x)),
        be = ~ mean(.x)
      )
    )) %>% 
  # Ungroup the data
  ungroup() %>% 
  # Centering the scaling group-level predictors
  mutate(
    across(
      ends_with("_be"),
      ~ (.x - mean(.x))/(2*sd(.x))
    ),
    # Declare country as a factor
    country = as_factor(country)
  )

# Read in the weights matrix
ipw_mat <- read_rds("data/ipw_weights_full.rds")

#------------------------------------------------------------------------------#
#---------------------Estimate the Naive Outcome Stage Model--------------------
#------------------------------------------------------------------------------#

# Naive Outcome Model
naive_outcome_model <- bf(
  y_treated_conf ~ treat_bin_conf + beta_wi + beta_be + gamma_wi + gamma_be + 
    delta_wi + delta_be + (1 | country) + ar(time = time, gr = country, p = 1),
  center = FALSE,
  family = gaussian(link = "identity")
)

# Priors for b parameters are based on (sd(y)/sd(x))*2.5
outcome_priors <- prior(normal(0, 18.57), class = "b", coef = "beta_wi") +
  prior(normal(0, 18.57), class = "b", coef = "beta_be") +
  prior(normal(0, 10.43), class = "b", coef = "gamma_wi") +
  prior(normal(0, 10.43), class = "b", coef = "gamma_be") +
  prior(normal(0, 18.38), class = "b", coef = "delta_wi") +
  prior(normal(0, 18.38), class = "b", coef = "delta_be") +
  prior(normal(0, 2.5), class = "b", coef = "treat_bin_conf") +
  prior(normal(0, 50), class = "b", coef = "Intercept") +
  prior(exponential(0.01518), class = "sd") +
  prior(normal(0, 0.5), class = "ar")

# Fit the model using brms----
naive_outcome_fit <- brm(
  naive_outcome_model,
  data = sim_data_df,
  prior = outcome_priors,
  chains = 6, 
  cores = 6L,
  iter = 8000,
  seed = 12345, 
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  refresh = 100,
  save_model = str_c(stan_dir, "Naive_Outcome_LMM.stan"),
  file = str_c(fits_dir, "Naive_Outcome_LMM4"),
)

#------------------------------------------------------------------------------#
#----------Estimate the Latent IPTW Single-Weighted Outcome Stage Model---------
#------------------------------------------------------------------------------#

# Take advantage of brms functionality because I'm lazy
mod_data <- make_standata(naive_outcome_model, data = sim_data_df)

# Transpose the weights matrix for more efficient Stan access
weights <- t(ipw_mat)

# Prepare the data for use with Stan
iptw_stan_data <- list(
  N = nrow(sim_data_df),
  K = 8,
  Y = mod_data$Y,
  X = mod_data$X,
  J = mod_data$N_1,
  jj = mod_data$J_1,
  K_AR = 1,
  J_lag = mod_data$J_lag,
  IPW_N = dim(weights)[1],
  IPW = weights
)

# Load the Stan Model File for the Single Weighted Model
SWLMM_AR1 <- str_c(stan_dir, "IPTW_Outcome_SWLMM_AR1.stan")

# Compile the Stan model
swlmm_ar1_mod <- cmdstan_model(
  SWLMM_AR1, 
  dir = fits_dir,
  force_recompile = TRUE
)

# Print the model code
str_split(swlmm_ar1_mod$code(), pattern = ";", simplify = T)

# Fit the Outcome-Stage Model; This takes a little while to get started 
# since it has to calculate the mean and sd for 3,060 24k row columns
swlmm_ar1_fit <- swlmm_ar1_mod$sample(
  data = iptw_stan_data,
  seed = 123456,
  refresh = 50,
  output_dir = fits_dir,
  sig_figs = 5,
  parallel_chains = 6,
  chains = 6,
  iter_warmup = 5000,
  iter_sampling = 5000,
  max_treedepth = 11
)

# Write the model object to an RDS file
swlmm_ar1_fit$save_object(file = str_c(fits_dir, "IPTW_Outcome_SWLMM_AR1.rds"))

