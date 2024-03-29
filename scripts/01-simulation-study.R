#--------------------------Bayesian MSM Simulation------------------------------
#-Author: A. Jordan Nafa------------------------------Created: August 28, 2022-#
#-R Version: 4.1.2----------------------------------Revised: September 4, 2022-#

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
  "brms",
  "cmdstanr",
  "posterior",
  "bayesplot",
  "furrr",
  "ggblend",
  "tidybayes",
  "patchwork",
  install = FALSE
)

# Load the functions for the simulations
.helpers <- map(
  .x = list.files(
    path = "functions/",
    pattern = ".*.R",
    full.names = TRUE
  ),
  ~ source(.x)
)

#------------------------------------------------------------------------------#
#------------------------------Data Pre-Processing------------------------------
#------------------------------------------------------------------------------#

# Set the rng seed
set.seed(1234567)

# Simulate 2,000 datasets of varying dimensions 
reference_df <- expand.grid(
  groups = c(25, 45, 65, 85, 100),
  periods = c(20, 50),
  gamma = c(0, -0.5),
  treat_conf = FALSE,
  id = 1:100
) %>%
  mutate(
    # Calculate the sample size
    N = groups*periods,
    # Condition labels
    cond = if_else(gamma == 0, "Z Exogenous", "Z Endogenous")
  ) %>%
  # Nest the data by columns
  nest(sim_pars = c(groups, periods, treat_conf, gamma)) %>%
  # Simulate the datasets
  mutate(sim_data = map(
    .x = sim_pars,
    ~ dgp_sim(
      .groups = .x$groups, 
      .periods = .x$periods, 
      .true_gamma = .x$gamma, 
      .treat_conf = .x$treat_conf
    )
  )) %>%
  # Unnest the data dimensions
  unnest(cols = sim_pars)

#------------------------------------------------------------------------------#
#-------------------------------MSM Simulations---------------------------------
#------------------------------------------------------------------------------#

# Compile the Stan model
msm_sim_mod <- cmdstan_model("models/stan/IPTW_Outcome_Simulation.stan")

# Fit 3 models in parallel
plan(tweak(multisession, workers = 4))

# Estimate the weights for each model
msm_estimates <- reference_df %>%
  mutate(
    # Calculate the design stage msm for each dataset
    sim_data = map(
      .x = sim_data,
      ~ sim_msm_bayes_design(.x)
    ),
    # Build the data lists for the outcome model
    stan_data = future_map(
      .x = sim_data,
      .options = furrr_options(
        scheduling = 1,
        seed = TRUE,
        prefix = "prefix"
      ),
      .f = ~ make_msm_data(.x, shape_prior = c(2, 5)),
      .progress = TRUE
    ),
    # Calculate the design stage msm for each dataset
    outcome_estimates = future_map(
      .x = stan_data,
      .options = furrr_options(
        scheduling = 1,
        seed = TRUE,
        prefix = "prefix"
      ),
      .f = ~ sim_msm_bayes_outcome(.x, msm_sim_mod),
      .progress = TRUE
    ))

# Write the data frames with weights to a file because this takes 
# a long ass time to run
write_rds(msm_estimates, "data/msm_sims.rds")

#------------------------------------------------------------------------------#
#-------------------------------ARDL Simulations--------------------------------
#------------------------------------------------------------------------------#

# Compile the Stan model
ardl_sim_mod <- cmdstan_model("models/stan/ARDL_Simulation.stan")

# Fit 3 models in parallel
plan(tweak(multisession, workers = 4))

# Estimate the naive ARDL for each model
ardl_estimates <- msm_estimates %>%
  mutate(
    # Build the data lists for the outcome model
    ardl_stan_data = future_map(
      .x = sim_data,
      .options = furrr_options(
        scheduling = 1,
        seed = TRUE,
        prefix = "prefix"
      ),
      .f = ~ make_ardl_data(
        .x, 
        ardl.form = Y ~ X + X_Lag + Y_Lag + Y_Lag_2 + Z + Z_Lag,
        prior.scale = 1.5
      ),
      .progress = TRUE
    ),
    # Calculate the design stage msm for each dataset
    ardl_estimates = future_map(
      .x = ardl_stan_data,
      .options = furrr_options(
        scheduling = 1,
        seed = TRUE,
        prefix = "prefix"
      ),
      .f = ~ sim_ardl_bayes(.x, ardl_sim_mod),
      .progress = TRUE
    ))

## Write the updated version to a file
write_rds(ardl_estimates, "data/final_sims.rds")

#------------------------------------------------------------------------------#
#------------------------Build the data frame of Results------------------------
#------------------------------------------------------------------------------#

# Create a new tibble with the information for the graphs
sim_results <- ardl_estimates %>%
  transmute(
    across(c(id:periods, gamma)),
    MSM_X_Lag = NA_real_,
    MSM_X = NA_real_,
    ARDL_X_Lag = NA_real_,
    ARDL_Y_Lag = NA_real_,
    ARDL_X = NA_real_,
    truth = 0
  )

# Pull the posterior median for each model
for (i in seq_along(ardl_estimates$sim_data)) {
  # MSM Estimate for X[t-1]
  sim_results[i, "MSM_X_Lag"] <- ardl_estimates$outcome_estimates[[i]] %>% 
    .[4, "mean"]
  # MSM Estimate for X[t-1]
  sim_results[i, "MSM_X"] <- ardl_estimates$outcome_estimates[[i]] %>% 
    .[3, "mean"]
  # ARDL Estimate for X[t-1]
  sim_results[i, "ARDL_X_Lag"] <- ardl_estimates$ardl_estimates[[i]] %>% 
    .[.$variable == "b[2]", "mean"]
  # ARDL Estimate for Y[t-1]
  sim_results[i, "ARDL_Y_Lag"] <- ardl_estimates$ardl_estimates[[i]] %>% 
    .[.$variable == "b[4]", "mean"]
  # ARDL Estimate for X
  sim_results[i, "ARDL_X"] <- ardl_estimates$ardl_estimates[[i]] %>% 
    .[.$variable == "b[1]", "mean"]
}

# Calculate the ARDL Result
sim_results <- sim_results %>% 
  mutate(ARDL_Estimate = ARDL_X_Lag + ARDL_Y_Lag * ARDL_X)

## Write the updated version to a file
write_rds(sim_results, "data/sim_results.rds")

#------------------------------------------------------------------------------#
#--------------------------Verify Parameter Recovery----------------------------
#------------------------------------------------------------------------------#

# Calculate the loss function
error_by_groups <- sim_results %>% 
  # Group by dimensions
  group_by(id, groups, cond) %>% 
  # Calculate the mean
  summarise(across(
    c(MSM_X_Lag, ARDL_Estimate), 
    ~ sqrt(mean((truth - .x)^2)),
    .names = "{.col}_rmse"
  )) %>% 
  # Pivot the result to long form
  pivot_longer(cols = ends_with("rmse")) %>% 
  ## Set the names for the facets
  mutate(name = if_else(name == "MSM_X_Lag_rmse", "MSM Bias", "ARDL Bias"))

# Plot RMSE by number of groups
bias_groups_plot <- ggplot(error_by_groups, aes(x = groups, y = value)) +
  # Facet the plot by condition
  facet_wrap(~ cond) +
  geom_point(
    aes(fill = name, shape = name), 
    size = 3,
    position = "jitter"
  ) * blend("multiply") +
  scale_fill_manual(values = c("#9400D3", "#00CD00")) +
  scale_shape_manual(values = c(22, 23)) +
  scale_x_continuous(breaks = seq(20, 100, 20)) +
  scale_y_continuous(limits = c(0, 0.10)) +
  labs(
    y = latex2exp::TeX(r'($\sqrt{E(\hat{\theta} - \theta)^{2}}$)'),
    x = "Number of Groups",
    fill = "Condition",
    shape = "Condition",
  ) +
  theme_bw(base_family = "serif", base_size = 30) +
  theme(legend.position = "top")

# Calculate the loss function
error_by_periods <- sim_results %>% 
  # Group by dimensions
  group_by(id, periods, cond) %>% 
  # Calculate the mean
  summarise(across(
    c(MSM_X_Lag, ARDL_Estimate), 
    ~ sqrt(mean((truth - .x)^2)),
    .names = "{.col}_rmse"
  )) %>% 
  # Pivot the result to long form
  pivot_longer(cols = ends_with("rmse")) %>% 
  ## Set the names for the facets
  mutate(name = if_else(name == "MSM_X_Lag_rmse", "MSM Bias", "ARDL Bias"))

# Plot RMSE by number of periods
bias_periods_plot <- ggplot(error_by_periods, aes(x = periods, y = value)) +
  # Facet the plot by condition
  facet_wrap(~ cond) +
  geom_point(
    aes(fill = name, shape = name), 
    size = 3,
    position = "jitter"
  ) * blend("multiply") +
  scale_fill_manual(values = c("#9400D3", "#00CD00")) +
  scale_shape_manual(values = c(22, 23)) +
  scale_y_continuous(limits = c(0, 0.10)) +
  scale_x_continuous(breaks = c(20, 50)) +
  labs(
    y = latex2exp::TeX(r'($\sqrt{E(\hat{\theta} - \theta)^{2}}$)'),
    x = "Number of Periods",
    fill = "Condition",
    shape = "Condition",
  ) +
  theme_bw(base_family = "serif", base_size = 30) +
  theme(legend.position = "top")

# Patch the two plots into one
rmse_plots <- bias_periods_plot/bias_groups_plot +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
    )

# Write the plot to a file
ggsave(
  "rmse_lagx_sim_results.jpeg",
  plot = rmse_plots,
  device = "jpeg",
  path = "figures",
  width = 16,
  height = 16,
  units = "in",
  dpi = "retina",
  type = "cairo"
)

# Data for plotting the distributions of the estimates
param_results <- sim_results %>% 
  # Pivot estimates to long form
  pivot_longer(
    cols = c(MSM_X_Lag:MSM_X, ARDL_X, ARDL_Estimate),
    names_to = "param",
    values_to = "estimate"
  ) %>% 
  # Make facet labels
  mutate(
    parameter = case_when(
      param == "MSM_X_Lag" ~ latex2exp::TeX(r'(Distribution of MSM Estimates: $X_{it-1}$)', output = "character"),
      param == "ARDL_Estimate" ~ latex2exp::TeX(r'(Distribution of ARDL Estimates: $X_{it-1}$)', output = "character"),
      param == "MSM_X" ~ latex2exp::TeX(r'(Distribution of MSM Estimates: $X_{it}$)', output = "character"),
      param == "ARDL_X" ~ latex2exp::TeX(r'(Distribution of ARDL Estimates: $X_{it}$)', output = "character")
    ),
    model = case_when(
      str_detect(param, "MSM") ~ "MSM",
      str_detect(param, "ARDL") ~ "ARDL",
    ),
    truth = case_when(
      str_detect(param, "X_Lag|Estimate") ~ 0,
      str_detect(param, "MSM_X|ARDL_X") ~ -0.1
    )
  )

# Distribution of the means for the lag of x
x_lag_estimates <- param_results %>% 
  filter(truth == 0) %>% 
  ggplot(., aes(x = estimate, y = cond)) +
  facet_wrap(~ parameter, labeller = label_parsed) +
  stat_slabinterval(
    aes(slab_alpha = stat(pdf), fill = model, shape = model),
    fill_type = "gradient",
    point_interval = "mean_qi",
    show.legend = FALSE
  ) +
  geom_vline(xintercept = 0, lty = "dashed") +
  scale_shape_manual(values = c(22, 23)) +
  scale_fill_manual(values = palettetown::pokepal(2, 4)) +
  labs(
    y = "",
    x = "Posterior Mean Estimates"
  ) +
  theme_bw(base_family = "serif", base_size = 30) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

# Distribution of the means for x
x_contemp_estimates <- param_results %>% 
  filter(truth == -0.1) %>% 
  ggplot(., aes(x = estimate, y = cond)) +
  facet_wrap(~ parameter, labeller = label_parsed) +
  stat_slabinterval(
    aes(slab_alpha = stat(pdf), fill = model, shape = model),
    fill_type = "gradient",
    point_interval = "mean_qi",
    show.legend = FALSE
  ) +
  geom_vline(xintercept = -0.1, lty = "dashed") +
  scale_fill_manual(values = palettetown::pokepal(2, 4)) +
  scale_shape_manual(values = c(22, 23)) +
  labs(
    y = "",
    x = "Posterior Mean Estimates"
  ) +
  theme_bw(base_family = "serif", base_size = 30) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

# Combine the two plots in to one
dist_plots <- x_lag_estimates/x_contemp_estimates +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

# Write the plot to a file
ggsave(
  "sim_results_distributions.jpeg",
  plot = dist_plots,
  device = "jpeg",
  path = "figures",
  width = 16,
  height = 16,
  units = "in",
  dpi = "retina",
  type = "cairo"
)
