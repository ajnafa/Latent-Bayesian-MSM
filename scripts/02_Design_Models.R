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
  "cmdstanr",
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
#  1.19472   3.60030  -0.07034  -0.09988 

# Get the true values for the treatment
(true_resp <- sim_data_ls$true_resp)

# Check distribution of the treatment
xtabs(~ country + post_treat, sim_data_ls$data)

# Extract the data
sim_data_df <- sim_data_ls$data %>%
  # Group the data by country
  group_by(country) %>%
  # Centering the predictors
  mutate(across(
    c(time, beta:delta),
    list(
      wi = ~ (.x - mean(.x))/(2*sd(.x)),
      be = ~ mean(.x)
    )
  )) %>% 
  # Ungroup the data
  ungroup()

#------------------------------------------------------------------------------#
#------------------------Estimate the Design Stage Model------------------------
#------------------------------------------------------------------------------#

# Prepare the data for use with Stan
hlogit_stan_data <- list(
  N = nrow(sim_data_df),
  K = 4,
  Y = sim_data_df$post_treat,
  X = select(sim_data_df, ends_with("wi")),
  J = max(sim_data_df$country),
  jj = sim_data_df$country
)

# Load the Stan Model File
hlogit_prop <- str_c(stan_dir, "Design_Stage_HLogit_1.stan")

# Compile the Stan model
hlogit_prop_mod <- cmdstan_model(
  hlogit_prop, 
  dir = fits_dir,
  force_recompile = TRUE
)

# Print the model code
str_split(hlogit_prop_mod$code(), pattern = ";", simplify = T)

# Fit the Design-Stage Model; Run Time is 36.0 seconds
hlogit_prop_fit <- hlogit_prop_mod$sample(
  data = hlogit_stan_data,
  seed = 123456,
  refresh = 50,
  output_dir = fits_dir,
  sig_figs = 5,
  parallel_chains = 6,
  chains = 6,
  iter_warmup = 4000,
  iter_sampling = 2000,
  max_treedepth = 11
)

# Write the model object to an RDS file
hlogit_prop_fit$save_object(file = str_c(fits_dir, "Design_Stage_HLogit_1.rds"))

# Leave One Out Cross Validation
hlogit_prop_fit_loo <- hlogit_prop_fit$loo(cores = 4L)

