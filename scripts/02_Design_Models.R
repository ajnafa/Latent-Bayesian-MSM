#--------------Inverse Probability Weights for the Gaussian Model---------------
#-Author: A. Jordan Nafa-----------------------------Created: January 27, 2022-#
#-R Version: 4.1.0-----------------------------------Revised: January 27, 2022-#

# Set Project Options----
options(
  digits = 4, # Significant figures output
  scipen = 999, # Disable scientific notation
  repos = getOption("repos")["CRAN"]
)

# Load the necessary libraries----
pacman::p_load(
  "tidyverse",
  "data.table",
  "dtplyr",
  "brms",
  "tidybayes"
)

# Load the helper functions
source("scripts/00_Helper_Functions.R")

# Specify the directory for stan models
stan_dir <- "models/stan/"

# Specify the directory for fitted model files
fits_dir <- "models/fits/design_stage/"

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
sim_data_df <- sim_data$data %>%
  # Group the data by country
  group_by(country) %>%
  # Centering and scaling the predictors
  mutate(
    across(
      c(time, beta:delta),
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

#------------------------------------------------------------------------------#
#------------------------Estimate the Design Stage Model------------------------
#------------------------------------------------------------------------------#

# Specify the formula for the propensity model
bf_hlogit_prop_mod_re <- bf(
  treat_bin_conf ~ beta_wi + beta_be + gamma_wi + gamma_be + delta_wi + 
    delta_be + (1 | time) + (1 | country),
  center = FALSE,
  family = bernoulli(link = "logit")
)

# Priors for the model parameters
prop_priors <- prior(normal(0, 1), class = "b") +
  prior(student_t(4, 0, 1), class = "b", coef = "Intercept") +
  prior(exponential(0.8), class = "sd")

# Fit the model using brms----
hlogit_prop_re_fit <- brm(
  bf_hlogit_prop_mod_re,
  data = sim_data_df,
  prior = prop_priors,
  chains = 6, 
  cores = 6L,
  iter = 8000,
  seed = 12345, 
  backend = "cmdstanr",
  save_pars = save_pars(all = TRUE),
  refresh = 100,
  save_model = "Design_Stage_HLogit",
  file = str_c(fits_dir, "Design_Stage_HLogit"),
  str_c(stan_dir, "Design_Stage_HLogit.stan")
)

# Add LOO and Bayes R2 to the Model for the Full Data
hlogit_prop_re_fit <- add_criterion(
  hlogit_prop_re_fit,
  criterion = c("loo", "loo_R2"),
  cores = 4,
  seed = 666
)

# Generate posterior expectations for each observation
pred_probs_chains <- posterior_epred(hlogit_prop_re_fit)

# Rows are posterior draws, columns are original rows in dataset
dim(pred_probs_chains) # 24000 x 3060 matrix

# Transpose the matrix so that columns are posterior draws
ipw_matrix <- t(pred_probs_chains) %>%
  # Coerce the matrix to a tibble for manipulation
  as_tibble(.name_repair = "universal") %>% 
  # Add the treatment column for calculating inverse probability weights
  mutate(treat_bin_conf = sim_data_df$treat_bin_conf) %>% 
  # Calculate the inverse probability weights
  mutate(across(
    starts_with("..."),
    ~ (treat_bin_conf / .x) + ((1 - treat_bin_conf) / (1 - .x))
  )) %>% 
  # Get rid of treatment column
  select(-treat_bin_conf)

# Write the weights to a file
write_rds(ipw_matrix, "data/ipw_weights_full.rds")

