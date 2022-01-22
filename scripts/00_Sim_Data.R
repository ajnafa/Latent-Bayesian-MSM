#------------------------Data for the Simulation Study--------------------------
#-Author: A. Jordan Nafa------------------------------Created: January 1, 2022-#
#-R Version: 4.1.0------------------------------------Revised: January 1, 2022-#

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
  "truncnorm",
  "scales",
  "ids"
)

#------------------------------------------------------------------------------#
#----------------Simulating Data for the Gaussian Mixed Models------------------
#------------------------------------------------------------------------------#

# Defining the function for simulating cross-sectional time series data----
sim_xts <- function(J, T, seed = 12345, ...) {
  # This function is adapted from the example in Andrew's blog post here:
  # https://www.andrewheiss.com/blog/2020/12/03/ipw-tscs-msm/

  # Set the RNG seed
  set.seed(seed)

  # Define the Dimensions for the Simulated Data
  Countries <- J # Number of Countries
  Periods <- T # Number of Time Periods
  N <- Countries * Periods # Number of Observations

  # Generate the data frame to loop across
  pars_grid <- expand_grid(
    jj = seq(1, Countries, 1),
    tt = seq(1, Periods, 1)
  )

  # Simulate the year-specific shocks shared across countries
  for (t in 1:Periods) {
    # Population Time Shock ~ Normal+(0.1, 0.1)
    pars_grid$tt_pop_shock[pars_grid$tt == t] <- rtruncnorm(1, 0.1, 0.1, a = 0)
    # GDP Time Shock ~ Normal(0, 0.03)
    pars_grid$tt_gdp_shock[pars_grid$tt == t] <- rnorm(1, 0, 0.03)
    # Democracy Time Shock ~ Normal(0, 0.02)
    pars_grid$tt_democ_shock[pars_grid$tt == t] <- rnorm(1, 0, 0.02)
    # Corruption Time Shock ~ Normal(0, 0.01)
    pars_grid$tt_corr_shock[pars_grid$tt == t] <- rnorm(1, 0, 0.01)
  }

  # Simulate the baseline values by country
  for (j in 1:Countries) {
    # Baseline Population ~ Normal(16.5, 3)
    pars_grid$jj_pop_base[pars_grid$jj == j] <- rnorm(1, 16.5, 3)
    # Population Growth ~ Normal(0.04, 0.007)
    pars_grid$jj_pop_grw[pars_grid$jj == j] <- rnorm(Periods, 0.04, 0.007)
    # Baseline GDP ~ Uniform(16, 20.5)
    pars_grid$jj_gdp_base[pars_grid$jj == j] <- runif(1, 16, 20.5)
    # GDP Growth ~ Normal(0.05, 0.05)
    pars_grid$jj_gdp_grw[pars_grid$jj == j] <- rnorm(Periods, 0.05, 0.05)
    # Baseline Democracy ~ Uniform(0.2, 0.6)
    pars_grid$jj_democ_base[pars_grid$jj == j] <- runif(1, 0.2, 0.6)
    # Democracy Growth ~ Normal(0.03, 0.01)
    pars_grid$jj_democ_grw[pars_grid$jj == j] <- rnorm(Periods, 0.03, 0.01)
    # Baseline Corruption ~ Uniform(0.1, 0.6)
    pars_grid$jj_corr_base[pars_grid$jj == j] <- runif(1, 0.1, 0.6)
    # Corruption Growth ~ Normal(0, 0.04)
    pars_grid$jj_corr_grw[pars_grid$jj == j] <- rnorm(Periods, 0, 0.04)
    # Baseline Vacation Days ~ Uniform(5, 12)
    pars_grid$jj_pvcat_base[pars_grid$jj == j] <- runif(1, 5, 12)
    # Vacation Days Growth ~ Normal(0.75, 0.75)
    pars_grid$jj_pvcat_grw[pars_grid$jj == j] <- rnorm(Periods, 0.75, 0.75)
    # Baseline Policy Intervention ~ Uniform(5, 12)
    pars_grid$jj_policy_base[pars_grid$jj == j] <- runif(1, 25, 75)
    # Policy Intervention Growth ~ Normal(1.5, 1)
    pars_grid$jj_policy_grw[pars_grid$jj == j] <- rnorm(1, 1.5, 1)
    # Baseline Public Happiness ~ Uniform(5, 12)
    pars_grid$jj_happin_base[pars_grid$jj == j] <- runif(1, 20, 45)
    # Baseline Public Happiness ~ Normal(0.4, 0.05)
    pars_grid$jj_happin_grw[pars_grid$jj == j] <- runif(1, 20, 45)
  }

  # Construct the simulated dataset
  sim_df <- pars_grid %>%
    mutate(
      # Country Identifier
      country = jj,
      # Temporal Identifier
      year = tt,
      # Log Population
      ln_pop = jj_pop_base + ((tt * 0.1) * jj_pop_grw),
      # Log GDP
      ln_gdp = jj_gdp_base + (0.4 * ln_pop) +
        tt_gdp_shock + ((tt * 0.1) * jj_gdp_grw),
      # GDP
      gdp = exp(ln_gdp),
      # Total Population
      population = exp(ln_pop),
      # GDP Per Capita
      gdp_pcap = gdp / population,
      # Log GDP Per Capita
      ln_gdp_pcap = log(gdp_pcap),
      # Democracy
      democracy = jj_democ_base + tt_democ_shock + ((tt * 0.1) * jj_democ_grw),
      # Rescale Democracy to range from 0-100
      democ_scaled = rescale(democracy, to = c(0, 1)) * 100,
      # Corruption
      corruption = jj_corr_base + tt_corr_shock + ((tt * 0.1) * jj_corr_grw),
      # Rescale Democracy to range from 0-100
      corr_scaled = rescale(corruption, to = c(0, 1)) * 100,
      # Number of Vacation Days
      vacation_days = jj_pvcat_base + ((tt * 0.1) * jj_pvcat_grw) + (0.00012 * gdp_pcap) +
        (0.09 * democ_scaled) + (-0.12 * corr_scaled),
      # Round Vacation Days to the nearest integer
      vacation_days = round(vacation_days, 0),
      # The Response Variable, the Assumed Causal Effect of Vacation Days is 1.7
      happin_vacat = jj_happin_base + ((tt * 0.1) * jj_happin_grw) + (0.00015 * gdp_pcap) +
        (0.11 * democ_scaled) + (-0.15 * corr_scaled) + (1.7 * vacation_days),
      # Treatment + Outcome, 6-hour workday; Generate a latent score for policy adoption
      policy_score = jj_policy_base + ((tt * 0.5) * jj_policy_grw) + (0.00012 * gdp_pcap) +
        (0.15 * democ_scaled) + (-0.29 * corr_scaled) + rnorm(N, 0, 10),
      # Rescale it to range from 0.05 to 0.60 and use it as a bernoulli probability
      policy_prob = case_when(
        tt >= 3 ~ rescale(policy_score, to = c(0.05, 0.60)),
        TRUE ~ 0 # Ensure that countries don't adopt the policy before period 3
      ),
      # Whether a country adopted the policy
      policy_adoption = rbinom(N, 1, policy_prob)
    ) %>%
    # Group the data by country
    group_by(country) %>%
    # Code the policy intervention 1 for every year after adoption
    mutate(policy = ifelse(rleid(policy_adoption) > 1, 1, 0)) %>%
    # Ungroup the data
    ungroup() %>%
    # Create the outcome variable for the policy, the causal effect is +7 points
    mutate(
      happin_policy = (rnorm(n(), 7, 3) * policy) + jj_happin_base +
        ((tt * 0.1) * jj_happin_grw) + (0.00015 * gdp_pcap) + (0.11 * democ_scaled) +
        (-0.15 * corr_scaled) + (0.7 * vacation_days) + rnorm(n(), 8, 2)
    ) %>%
    # Regroup the data by country
    group_by(country) %>%
    # Generate lagged variables
    mutate(across(
      c(happin_policy, policy, vacation_days, happin_vacat),
      ~ lag(.x, 1L),
      .names = "lag_{.col}"
    )) %>%
    # Ungroup the data
    ungroup() %>%
    # Remove first year after making lags
    filter(year > 1) %>% 
    # Keep only the generated variables
    select(country:lag_happin_vacat)

  # Create a list with the compiled data sets
  sims_list <- list(
    "Simulated Data" = sim_df,
    "Simulated Parameters" = pars_grid
  )

  # Return just the created list object
  return(sims_list)
}

# Arguments: 
#  J       The number of countries to simulate data for.
#
#  T       The number of periods to simulate data for per country.
#
#  seed    Seed for the random number generator. Defaults to 12345
#
#  ...     Additional arguments for future development. Currently 
#          unused

# Simulate the data with 50 countries and 30 years
happiness_ls <- sim_xts(J = 50, T = 30, seed = 12345)