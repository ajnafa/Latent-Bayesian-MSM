#------------------------Data for the Simulation Study--------------------------
#-Author: A. Jordan Nafa------------------------------Created: January 1, 2022-#
#-R Version: 4.1.0------------------------------------Revised: January 1, 2022-#

# Set Project Options----
options(
  digits = 4, # Significant figures output
  scipen = 999, # Disable scientific notation
  repos = getOption("repos")["CRAN"],
  mc.cores = 12L
)

# Load the necessary libraries----
pacman::p_load(
  "tidyverse",
  "data.table",
  "dtplyr",
  "truncnorm",
  "scales"
)

# Load the helper functions
source("scripts/00_Helper_Functions.R")

#------------------------------------------------------------------------------#
#----------------Simulating Data for the Gaussian Mixed Models------------------
#------------------------------------------------------------------------------#

# Defining the function for simulating cross-sectional time series data----
sim_AR1_xts <- function(J, T, resp_coefs = NULL, 
                        treat_coefs = NULL, seed = 12345, conf, ...) {
  # Define the Dimensions for the Simulated Data
  Countries <- J # Number of Countries
  Periods <- T # Number of Time Periods
  N <- Countries * Periods # Number of Observations
  
  # Set the RNG seed
  set.seed(seed)
  
  # Check if user-specified covariates for the response are provided----
  if (is.null(resp_coefs)) {
    # If no coefficient values are provided simulate from a normal dist
    resp_coef_vals <- sample(rnorm(N, 0, 3), 5)
    names(resp_coef_vals) <- c("intercept", "treat", "beta", "gamma", "delta")
  }
  # If user covariate values of length 5 are provided use those
  else if (!is.null(resp_coefs) & length(resp_coefs) == 5) {
    resp_coef_vals <- resp_coefs
    names(resp_coef_vals) <- c("intercept", "treat", "beta", "gamma", "delta")
  }
  # Throw an error is length of provided vector is wrong
  else {
    stop("User-specified covariates for the response must be NULL 
         or a numeric vector of length 4")
  }
  
  # Check if user-specified covariates for the treatment are provided----
  if (is.null(treat_coefs)) {
    # If no coefficient values are provided simulate from a normal dist
    treat_coef_vals <- sample(rnorm(N, 0, 1), 4)
    names(treat_coef_vals) <- c("intercept", "beta", "gamma", "delta")
  }
  # If user covariate values of length 4 are provided use those
  else if (!is.null(treat_coefs) & length(treat_coefs) == 4) {
    treat_coef_vals <- treat_coefs
    names(treat_coef_vals) <- c("intercept", "beta", "gamma", "delta")
  }
  # Throw an error is length of provided vector is wrong
  else {
    stop("User-specified covariates for the response must be NULL 
         or a numeric vector of length 4")
  }
  
  # Generate the data frame to loop across
  resp_grid <- expand_grid(
    country = seq(1, Countries, 1),
    time = seq(1, Periods, 1)
  )
  
  # For y[1[j]] ~ half-normal(5, 10)
  resp_grid$Y[resp_grid$time == 1] <- 
    rtruncnorm(Countries, a = 0, mean = 5, sd = 10)
  
  # Transpose the data into a j x t matrix
  resp_mat <- resp_grid %>% 
    # Pivot the data to wide form
    pivot_wider(
      id_cols = time,
      values_from = Y,
      names_from = country
    ) %>% 
    # Coerce the tibble to a matrix
    as.matrix()
  
  # Loop across countries and simulate the baseline value for the AR(1) process----
  for (j in 2:(Countries + 1)) {
    # For each y[2:T, j] ~ half-normal(mu[t - 1, j], 3)
    for (t in 2:Periods) {
      resp_mat[t, j] <- 
        rtruncnorm(1, a = 0, mean = resp_mat[t - 1, j], sd = 3)
    }
  }
  
  # Coerce the matrix back to a data frame----
  resp_grid <- as_tibble(resp_mat) %>% 
    # Convert the data back to long form
    pivot_longer(
      cols = 2:dim(resp_mat)[2],
      values_to = "Y",
      names_to = "country"
    ) %>% 
    # Sort the data by country and time
    arrange(country, time) %>% 
    # Coerce identifiers back to numeric
    mutate(across(time:country, ~ as.numeric(.x)))
  
  # Simulate country-specific time invariant confounding
  for (j in 1:Countries) {
    resp_grid$Conf[resp_grid$country == j] <- runif(1, -conf, conf)
  }
  
  # Simulate the data for the fixed effects parameters----
  fixeffs_grid <- expand_grid(
    country = seq(1, Countries, 1),
    time = seq(1, Periods, 1)
  ) %>% 
    # Coerce the tibble to a data table for easier indexing
    as.data.table()
  
  # beta[1[j]] ~ Normal(5, 2.5) + Noise
  fixeffs_grid$beta[fixeffs_grid$time == 1] <- 
    rnorm(Countries, 8, 2.5) + rnorm(Countries, 0, 1)
  # gamma[1[j]] ~ Uniform(16.5, 20.5) + Noise
  fixeffs_grid$gamma[fixeffs_grid$time == 1] <- 
    runif(Countries, 16.5, 20.5) + rnorm(Countries, 0, 1)
  # delta[1[j]] ~ Normal(10, 1) + Noise
  fixeffs_grid$delta[fixeffs_grid$time == 1] <- 
    rnorm(Countries, 10, 1) + rnorm(Countries, 0, 1)
  
  # Simulate the evolution of the data for the coefficients over time
  for (j in 1:Countries) {
    for (t in 2:Periods) {
      # beta[2:T[j]] ~ Normal(beta[t - 1], 1)
      fixeffs_grid[time == t & country == j]$beta <- 
        rnorm(1, fixeffs_grid[time == (t - 1) & country == j]$beta, sd = 1)
      # gamma[2:T[j]] ~ Normal(gamma[t - 1], 2.5)
      fixeffs_grid[time == t & country == j]$gamma <- 
        rnorm(1, fixeffs_grid[time == (t - 1) & country == j]$gamma, sd = 2.5)
      # delta[2:T[j]] ~ Normal(delta[t - 1], 0.25)
      fixeffs_grid[time == t & country == j]$delta <- 
        rnorm(1, fixeffs_grid[time == (t - 1) & country == j]$delta, sd = 1.25)
    }
  }

  # Construct the full data
  treatment_grid <- resp_grid %>%
    # Combine the two data frames
    left_join(fixeffs_grid, by = c("time", "country")) %>%
    # Create the treatment process
    mutate(
      # Model the treatment assignment as a function of time and the covariates
      treat = treat_coef_vals["intercept"] + (0.0005*(time - 2)) +
        (treat_coef_vals["beta"] * beta) + (treat_coef_vals["gamma"] * gamma) +
        (treat_coef_vals["delta"] * delta) + rnorm(N, 0, 0.5),
      # Generate confounded treatment
      treat_conf = treat + Conf,
      # Rescale the latent treatment
      treat_scaled = rescale(treat, to = c(0, 1)),
      # Rescale the confounded latent treatment
      treat_scaled_conf = rescale(treat_conf, to = c(0, 1)),
      # Generate a binary event indicator
      treat_bin = rbinom(N, 1, plogis(treat_scaled)),
      # Generate a confounded binary event indicator
      treat_bin_conf = rbinom(N, 1, plogis(treat_scaled_conf)),
      # Recode Time
      time = time - 1
    ) %>%
    # Group the data by country
    group_by(country) %>%
    # Induce dependence of the response on the treatment
    mutate(
      # Confounded Response and Unconfounded Treatment
      y_conf_treated = resp_coef_vals["intercept"] + (resp_coef_vals["beta"] * beta) +
        (resp_coef_vals["gamma"] * gamma) + (resp_coef_vals["delta"] * delta) +
        (resp_coef_vals["treat"] * treat_bin) + rnorm(n(), 0, 1) + Conf,
      # Unconfounded Response and Confounded Treatment
      y_treated_conf = resp_coef_vals["intercept"] + (resp_coef_vals["beta"] * beta) +
        (resp_coef_vals["gamma"] * gamma) + (resp_coef_vals["delta"] * delta) +
        (resp_coef_vals["treat"] * treat_bin_conf) + rnorm(n(), 0, 1),
      # Confounded Response and Unconfounded Treatment
      y_conf_treated_conf = resp_coef_vals["intercept"] + (resp_coef_vals["beta"] * beta) +
        (resp_coef_vals["gamma"] * gamma) + (resp_coef_vals["delta"] * delta) +
        (resp_coef_vals["treat"] * treat_bin_conf) + rnorm(n(), 0, 1) + Conf
    ) %>%
    # Ungroup the data
    ungroup() %>%
    # Drop missing values induced by lags
    drop_na()
  
  # Create a list containing the data and true values
  sim_list <- list(
    "data" = treatment_grid,
    "true_resp" = resp_coef_vals,
    "true_treat" = treat_coef_vals
  )

  # Return the list object
  return(sim_list)
}

# Arguments: 
#  J             The number of countries to simulate data for.
#
#  T             The number of periods to simulate data for per country.
#
#  seed          Seed for the random number generator. Defaults to 12345
#
#  resp_coefs    Optional user supplied true values for the response
#                coefficients. If supplied, input must be a numeric vector
#                of length 4 of the form `c(treat, beta, gamma, delta)`.
#                See description below for more information.
#
#  treat_coefs   Optional user supplied true values for the treatment process
#                coefficients. If supplied, input must be a numeric vector
#                of length 4 of the form `c(intercept, beta, gamma, delta)`.
#                See description below for more information.
#
#  conf          A numeric value for the magnitude of time invariant 
#                confounding. Value is passed to `runif(-conf, conf)`
#
#  ...           Additional arguments for future development. Currently 
#                unused

# Description:
# The `sim_AR1_xts` function generates a simulated panel with J groups
# observed at T equally spaced periods. The assumed base data generation 
# process for the treatment assignment is:
#   X_{4tj} = \alpha + \0.001 \cdot time_{tj} + Z_{1tj}\beta + 
#       Z_{2tj}\gamma + Z_{3tj}\delta + \upsilon
# where \upsilon \sim N(0, \sigma_{j})
#
# The baseline values for the response vector are simulated based on
# an AR(1) correlation structure. The observed response is a function
# of the stochastic evolution of the AR(1) process, the intervention,
# and a three time varying variables, beta, delta, and gamma. The
# assumed causal effect of the treatment on the response is
#   theta_{tj} = \alpha + y_{t-1,j}\mu + X_{1tj}\beta + X_{2tj}\delta + 
#         X_{3tj}\gamma + X_{4tj}\nu + \upsilon + xi
# where \upsilon \sim N(0, \sigma_{j}) and \xi \sim N(0, 1)
#
# If no user-supplied values are specified for the response coefficients
# values are randomly drawn from a normal distribution with mean 0 and
# sd 3. Likewise, if no user-supplied values are specified for the treatment
# process coefficients, values are randomly drawn from a standard normal 
# distribution.
#
# The function returns a list object containing the complete simulated data 
# frame, the "true" values for the response, and the "true" values for the
# treatment process

# Examples: 
# Simulate the data with 50 countries and 100 years and randomly
# generated coefficient values for the response and treatment
# sim_data <- sim_AR1_xts(J = 50, T = 100, seed = 12345)

# Set the rng seed
set.seed(12345)

# Generate a vector of true values for the response coefficients
treat_trues <- rnorm(4, 0, 1)*rgamma(4, 1, 1)

# Simulate the data with 30 countries and 102 years and user-specified
# coefficient values for the response and treatment. I use 102 here
# so that the remaining number of periods after creating the lags is 100
sim_data <- sim_AR1_xts(
  J = 30, 
  T = 102, 
  seed = 12345,
  resp_coefs = c(12.3591, 6.8563, 3.1879, -2.4487, 5.4057),
  treat_coefs = treat_trues,
  conf = 15
  )

# Check distribution of the unconfounded treatment
xtabs(~ country + treat_bin, sim_data$data)

# Check distribution of the confounded treatment
xtabs(~ country + treat_bin_conf, sim_data$data)

# Print the true values of the coefficients in the response
sim_data$true_resp

# Plot a subset of the simulated data
sim_data$data %>% 
  # Extract a subset of the data
  filter(country %in% 1:9) %>% 
  # Drop missing values in time
  drop_na() %>% 
  # Initiate a ggplot object
  ggplot(aes(x = time, y = y_conf_treated, color = as_factor(treat_bin))) +
  # Facet the plot by group
  facet_wrap(~ country, scales = "free_y") + 
  # Add points for each value of y
  geom_point(size = 1.25) +
  # Customize the colors
  scale_color_manual(values = c("firebrick", "royalblue")) +
  # Set plot labels
  labs(
    x = "Time",
    y = "Response",
    title = "Simulated AR(1) Panel Data with Confounded Response and Unconfounded Treatment",
    color = "Treated"
  ) +
  # Set plot theme
  plot_theme(plot.margin = margin(3, 5, 5, 5, "mm"))

# Plot a subset of the simulated data with confounding
sim_data$data %>% 
  # Extract a subset of the data
  filter(country %in% 1:9) %>% 
  # Drop missing values in time
  drop_na() %>% 
  # Initiate a ggplot object
  ggplot(aes(x = time, y = y_treated_conf, color = as_factor(treat_bin_conf))) +
  # Facet the plot by group
  facet_wrap(~ country, scales = "free_y") + 
  # Add points for each value of y
  geom_point(size = 1.25) +
  # Customize the colors
  scale_color_manual(values = c("firebrick", "royalblue")) +
  # Set plot labels
  labs(
    x = "Time",
    y = "Response",
    title = "Simulated AR(1) Panel Data with Unconfounded Response and Confounded Treatment",
    color = "Treated"
  ) +
  # Set plot theme
  plot_theme(plot.margin = margin(3, 5, 5, 5, "mm"))

# Plot a subset of the simulated data with confounding
sim_data$data %>% 
  # Extract a subset of the data
  filter(country %in% 1:9) %>% 
  # Drop missing values in time
  drop_na() %>% 
  # Initiate a ggplot object
  ggplot(aes(x = time, y = y_conf_treated_conf, color = as_factor(treat_bin_conf))) +
  # Facet the plot by group
  facet_wrap(~ country, scales = "free_y") + 
  # Add points for each value of y
  geom_point(size = 1.25) +
  # Customize the colors
  scale_color_manual(values = c("firebrick", "royalblue")) +
  # Set plot labels
  labs(
    x = "Time",
    y = "Response",
    title = "Simulated AR(1) Panel Data with Confounded Response and Confounded Treatment",
    color = "Treated"
  ) +
  # Set plot theme
  plot_theme(plot.margin = margin(3, 5, 5, 5, "mm"))
  
# Write the simulated data to a file
write_rds(sim_data, "data/gaussian_ar1_sims.rds")
