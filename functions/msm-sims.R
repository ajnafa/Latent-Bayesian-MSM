# Function for the ADL models----
sim_ardl_bayes <- function(sim.data, ...) {
  
  # Define the model formula
  bayes_adl_form <- bf(
    Y ~ X + X_Lag + Y_Lag + Y_Lag_2 + Z + Z_Lag,
    family = brmsfamily(family = "gaussian", link = "identity")
  )
  
  # Priors for the model
  bayes_adl_priors <- prior(normal(0, 1.5), class = "b") +
    prior(normal(mean(Y), 1.5*sd(Y)), class = "Intercept") +
    prior(exponential(1/sd(Y)), class = "sigma")
  
  # ARDL (t-1) using brms
  ardl_bayes <- brm(
    formula = bayes_adl_form,
    data = sim.data,
    prior = bayes_adl_priors,
    chains = 4, 
    cores = 4L,
    iter = 6000,
    seed = 12345, 
    backend = "cmdstanr",
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.9),
    refresh = 0
  )
  
  # Calculate a summary of the draws
  ardl_result <- summarise_draws(ardl_bayes, .cores = 6L)
  
  # Return the data frame of draws
  return(ardl_result)
}

# Function for calculating the stabilized weights----
bayes_ipwt <- function(sim.data, psnum, psdenom) {
  # Generate posterior expectations for the numerator model
  preds_num <- t(brms::posterior_epred(psnum))
  
  # Generate posterior expectations denominator model
  preds_denom <- t(brms::posterior_epred(psdenom))
  
  # Calculate the numerator of the stabilized weights
  num_scores <- preds_num * sim.data$X + (1 - preds_num) * (1 - sim.data$X)
  
  # Calculate the denominator of the stabilized weights
  denom_scores <- preds_denom * sim.data$X + (1 - preds_denom) * (1 - sim.data$X)
  
  # Calculate the weights
  wts <- num_scores/denom_scores
  
  ## Coerce the output to a tibble
  weights <- tibble::as_tibble(wts, .name_repair = "minimal")
  
  ## Assign column names to each draw
  colnames(weights) <- paste("draw_", 1:ncol(weights), sep = "")
  
  ## TODO: Figure out the data.table code for this
  weights <- weights |>
    ## Add group and time columns
    dplyr::mutate(
      unit = sim.data$unit, 
      time = sim.data$time,
      .before = 1
    ) |>
    ## Group the data by identifier
    dplyr::group_by(unit) |>
    ## Calculate the cumulative product of the weights by unit
    dplyr::mutate(dplyr::across(
      dplyr::starts_with("draw"),
      ~ cumprod(tidyr::replace_na(.x, 1)))
    ) |>
    ## Ungroup the Data
    dplyr::ungroup()
  
  # Generate the Lags
  sim.data[, `:=` (
    wts_mean = rowMeans(wts),
    wts_sd = apply(wts, 1, sd),
    cws_mean = rowMeans(weights[, 3:ncol(weights)]),
    cws_med = apply(weights[, 3:ncol(weights)], 1, median),
    cws_sd = apply(weights[, 3:ncol(weights)], 1, sd),
    num_prob = rowMeans(preds_num),
    denom_prob = rowMeans(preds_denom))
  ]
  
  # Return the data frame with the weights info
  return(sim.data)
}

# Function for the design stage of the MSM models----
sim_msm_bayes_design <- function(sim.data, ...) {
  
  # Define the model formula for the design stage numerator
  bayes_num_form <- bf(
    X ~ X_Lag,
    family = brmsfamily(family = "bernoulli", link = "logit")
  )
  
  # Priors for the design stage model numerator
  bayes_num_priors <- prior(normal(0, 1), class = "b") +
    prior(normal(0, 2), class = "Intercept")
  
  # Design stage model numerator
  ps_num_bayes <- brm(
    formula = bayes_num_form,
    data = sim.data,
    prior = bayes_num_priors,
    chains = 4, 
    cores = 4L,
    iter = 6000,
    backend = "cmdstanr",
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.9),
    refresh = 0,
    ...
  )
  
  # Define the model formula for the design stage denominator
  bayes_denom_form <- bf(
    X ~ Y_Lag + Z + X_Lag,
    family = brmsfamily(family = "bernoulli", link = "logit")
  )
  
  # Priors for the design stage model denominator
  bayes_denom_priors <- prior(normal(0, 1), class = "b") +
    prior(normal(0, 2), class = "Intercept")
  
  # Design stage model denominator
  ps_denom_bayes <- brm(
    formula = bayes_denom_form,
    data = sim.data,
    prior = bayes_denom_priors,
    chains = 4, 
    cores = 4L,
    iter = 6000,
    backend = "cmdstanr",
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.9),
    refresh = 0,
    ...
  )
  
  # Calculating the stabilized weights
  out <- bayes_ipwt(
    sim.data, 
    psnum =  ps_num_bayes, 
    psdenom = ps_denom_bayes
  )
  
  # Return the updated data
  return(out)
}

# Function for the outcome stage of the MSM models----
sim_msm_bayes_outcome <- function(sim.data, .seed, ...) {
  
  # Define the model formula for the design stage numerator
  bayes_msm_form <- bf(Y ~ X + X_Lag, family = gaussian())
  
  # Take advantage of brms functionality because I'm lazy
  msm_data <- make_standata(bayes_msm_form, data = sim.data)
  
  # Prepare the data for use with Stan
  msm_data <- list(
    N = nrow(sim.data),
    K = msm_data$K,
    Y = msm_data$Y,
    X = msm_data$X,
    ipw_mu = sim.data$cws_mean,
    ipw_sigma = sim.data$cws_sd
  )
  
  # Compile the Stan model
  msm_sim_mod <- cmdstan_model(
    "models/stan/IPTW_Outcome_Simulation.stan",
    force_recompile = TRUE
  )
  
  # Fit the Outcome-Stage Model; Run time is approximately 10 minutes
  msm_sim_fit <- msm_sim_mod$sample(
    data = msm_data,
    seed = 123456,
    refresh = 0,
    sig_figs = 5,
    parallel_chains = 4,
    chains = 4,
    iter_warmup = 4000,
    iter_sampling = 4000,
    max_treedepth = 12,
    adapt_delta = 0.9
  )
  
  # Calculate a summary of the draws
  msm_result <- summarise_draws(msm_sim_fit, .cores = 6L)
  
  # Return the update data frame
  return(msm_result)
}

