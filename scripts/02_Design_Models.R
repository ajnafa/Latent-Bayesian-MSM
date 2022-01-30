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
  "tidybayes",
  "latex2exp"
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
sim_data_df <- sim_data_ls$data %>%
  # Group the data by country
  group_by(country) %>%
  # Centering and scaling the predictors
  mutate(
    across(
      beta:delta,
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
prop_priors <- prior(normal(0, 1.5), class = "b") +
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
  # Add Variable needed to calculating and structuring the weights
  mutate(
    # Time Identifier
    time = hlogit_prop_re_fit$data$time,
    # Country Identifier
    country = hlogit_prop_re_fit$data$country,
    # Add the treatment column for calculating inverse probability weights
    treat_bin_conf = hlogit_prop_re_fit$data$treat_bin_conf,
    # Observation ID
    obs = 1:n()
    ) %>% 
  # Calculate the inverse probability weights
  mutate(across(
    starts_with("..."),
    ~ (treat_bin_conf / .x) + ((1 - treat_bin_conf) / (1 - .x)),
    .names = "iptw_{.col}"
  )) %>% 
  # Pivot the data to wide form
  pivot_longer(cols = starts_with(c("...", "iptw"))) %>% 
  # Recode value type names
  mutate(name = case_when(
    str_detect(name, "iptw") ~ "iptw",
    TRUE ~ "prob"
  )) %>% 
  # group the data by identifiers
  group_by(across(time:name)) %>% 
  # Calculate the means and sds for each observation weight
  summarise(mu = mean(value), sigma = sd(value)) %>% 
  # Ungroup the data
  ungroup() %>% 
  # Pivot the collapsed data to wide form
  pivot_wider(
    id_cols = time:obs, 
    values_from = mu:sigma, 
    names_from = name
    )

# Write the weights to a file
write_rds(ipw_matrix, "data/iptw_weights_descs.rds")

#------------------------------------------------------------------------------#
#-----------------Plot of the Treated and Untreated Populations-----------------
#------------------------------------------------------------------------------#

# Plot of the propensity for treatment distributions
observed_pops <- ggplot() + 
  # Add a density slab for the observed probability of treated units
  stat_slab(
    data = filter(ipw_matrix, treat_bin_conf == 1), 
    aes(x = mu_prob, slab_alpha = stat(2*pdf)),
    fill = "#389078",
    side = "top"
    ) +
  # Add a density slab for the observed probability of untreated units
  stat_slab(
    data = filter(ipw_matrix, treat_bin_conf == 0), 
    aes(x = mu_prob, slab_alpha = stat(2*pdf)),
    fill = "#B82820",
    side = "bottom"
  ) +
  # Label for the treated slab
  annotate(
    geom = "label", 
    x = 0.52, 
    y = 0.5, 
    label = "Treated", 
    fill = "#389078", 
    color = "white", 
    hjust = 0,
    family = "serif",
    size = 6
    ) +
  # Label for the untreated slab
  annotate(
    geom = "label", 
    x = 0.52, 
    y = -0.5, 
    label = "Untreated", 
    fill = "#B82820", 
    color = "white", 
    hjust = 0,
    family = "serif",
    size = 6
  ) +
  # Adjust the number of breaks on the x axis
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  # Add custom theme settings
  plot_theme(plot.margin = margin(5, 1, 3, 4, "mm")) +
  # Add labels to the plot
  labs(
    x = TeX(r'($E\[Y_{tj}  \, | \, X_{n}\beta_{k} + \upsilon_{j} + \delta_{t}\]$)'),
    y = "Density",
    title = "Conditional Probability of Observed Treatment"
  ) +
  # Setting the parameters for the plot legend
  guides(slab_alpha = "none")

# Save the generated plot object as a .jpeg file
ggsave(
  filename = "Figure_1A_Population_Balance.jpeg",
  plot = observed_pops,
  device = "jpeg",
  path = "figures/",
  width = 16,
  height = 10,
  units = "in",
  dpi = "retina",
  type = "cairo",
  limitsize = FALSE
)

# Plot observed and pseudo-population distriubtions together
observed_pseudopops <- ggplot() + 
  # Treated Distribution (Pseudo-Population)
  geom_histogram(
    data = filter(ipw_matrix, treat_bin_conf == 1), 
    aes(x = mu_prob, weight = mu_iptw),
    bins = 50,
    fill = palettetown::pokepal(1)[4],
    color = palettetown::pokepal(1)[5]
    ) + 
  # Untreated Distribution (Pseudo-Population)
  geom_histogram(
    data = filter(ipw_matrix, treat_bin_conf == 0), 
    aes(x = mu_prob, weight = mu_iptw, y = -..count..), 
    bins = 50,
    fill = palettetown::pokepal(1)[12],
    color = palettetown::pokepal(1)[13]
  ) +
  # Treated Distribution (Observed)
  geom_histogram(
    data = filter(ipw_matrix, treat_bin_conf == 1), 
    aes(x = mu_prob), 
    bins = 50,
    fill = palettetown::pokepal(1)[1],
    color = palettetown::pokepal(1)[2]
  ) + 
  # Untreated Distribution (Observed)
  geom_histogram(
    data = filter(ipw_matrix, treat_bin_conf == 0), 
    aes(x = mu_prob, y = -..count..), 
    bins = 50,
    fill = palettetown::pokepal(1)[9],
    color = palettetown::pokepal(12)[11]
  ) +
  # Label for the treated observations
  annotate(
    geom = "label", 
    x = 0.45, 
    y = 65, 
    label = "Treated (Observed)", 
    fill = "#389078", 
    color = "white", 
    hjust = 0,
    family = "serif",
    size = 6
  ) +
  # Label for the treated pseudo-population
  annotate(
    geom = "label", 
    x = 0.45, 
    y = 90, 
    label = "Treated (Pseudo-Population)", 
    fill = "#98D048", 
    color = "white", 
    hjust = 0,
    family = "serif",
    size = 6
  ) +
  # Label for the untreated observations
  annotate(
    geom = "label", 
    x = 0.45, 
    y = -65, 
    label = "Untreated (Observed)", 
    fill = "#B82820", 
    color = "white", 
    hjust = 0,
    family = "serif",
    size = 6
  ) +
  # Label for the untreated pseudo-population
  annotate(
    geom = "label", 
    x = 0.45, 
    y = -90, 
    label = "Untreated (Pseudo-Population)", 
    fill = "#F86860", 
    color = "white", 
    hjust = 0,
    family = "serif",
    size = 6
  ) +
  # Adjust the number of breaks on the x axis
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  # Add custom theme settings
  plot_theme(plot.margin = margin(5, 1, 3, 4, "mm")) +
  # Add labels to the plot
  labs(
    x = TeX(r'($E\[Y_{tj} \, | \, X_{n}\beta_{k} + \upsilon_{j} + \delta_{t}\]$)'),
    y = "",
    title = "Observed Treatment and Pseudo-Population Distributions"
  ) +
  # Setting the parameters for the plot legend
  guides(slab_alpha = "none")
  
# Save the generated plot object as a .jpeg file
ggsave(
  filename = "Figure_1A_Pseudo-Population_Balance.jpeg",
  plot = observed_pseudopops,
  device = "jpeg",
  path = "figures/",
  width = 16,
  height = 10,
  units = "in",
  dpi = "retina",
  type = "cairo",
  limitsize = FALSE
)
