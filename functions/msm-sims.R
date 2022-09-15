#' Function for Calculating the Stabilized Inverse Probability Weights
#'
#' @param sim.data The simulated data table returned by `dgp_sim`
#' 
#' @param psnum The numerator model of class `brmsfit` for the stabilized 
#' weights
#' 
#' @param psdenom The denominator model of class `brmsfit` for the stabilized 
#' weights
#'
#' @return Returns a data table with the location and scale of the stabilized
#' inverse probability of treatment weights for the ATE of a binary treatment
#' 
#' @export bayes_ipwt
#'
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

# Function for Fitting the Design Stage Models of the MSM
#'
#' @param sim.data The simulated data table returned by `dgp_sim`
#' 
#' @param ... Additional arguments passed to `brms::brm`
#'
#' @return The original data table with the location and scale of the stabilized
#' inverse probability of treatment weights for the ATE of a binary treatment
#' 
#' @export sim_msm_bayes_design
#' 
sim_msm_bayes_design <- function(sim.data, ...) {

  # Design stage model numerator
  ps_num_bayes <- brm(
    formula = bf(X ~ X_Lag),
    data = sim.data,
    family = bernoulli(link = "logit"),
    prior = prior(normal(0, 1), class = "b") +
      prior(normal(0, 2), class = "Intercept"),
    chains = 4, 
    cores = 4L,
    iter = 4000,
    backend = "cmdstanr",
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.9, max_treedepth = 12),
    refresh = 0,
    ...
  )
  
  # Design stage model denominator
  ps_denom_bayes <- brm(
    formula = bf(X ~ Y_Lag + Z + X_Lag),
    data = sim.data,
    family = bernoulli(link = "logit"),
    prior = prior(normal(0, 1), class = "b") +
      prior(normal(0, 2), class = "Intercept"),
    chains = 4, 
    cores = 4L,
    iter = 4000,
    backend = "cmdstanr",
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = 0.9, max_treedepth = 12),
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

#' Function for Building the Data to Pass to the MSM Stan Model
#'
#' @param sim.data A data table object with the stabilized weight vectors 
#' returned by the `sim_msm_bayes_design` function
#' 
#' @param shape_prior A numeric vector of length 2 containing the location and
#' scale to use for the beta prior on the scale of the weights
#' 
#' @param ... Currently unused
#'
#' @return Returns a list object with the data to be passed to the Stan model
#' 
#' @export make_msm_data
#'
make_msm_data <- function(sim.data, shape_prior, ...) {
  
  # Take advantage of brms functionality because I'm lazy
  msm_data <- brms::make_standata(
    Y ~ X + X_Lag, 
    family = gaussian(), 
    data = sim.data
  )
  
  # Prepare the data for use with Stan
  msm_data <- list(
    N = nrow(sim.data),
    K = msm_data$K,
    Y = msm_data$Y,
    X = msm_data$X,
    ipw_mu = sim.data$cws_mean,
    ipw_sigma = sim.data$cws_sd,
    sd_prior_shape1 = shape_prior[1],
    sd_prior_shape2 = shape_prior[2]
  )
  
  return(msm_data)
}

#' Function for Fitting the MSM Outcome Model Simulations
#' 
#' @param msm.data A list object containing the data to be passed to the
#' Stan model as returned by the `make_msm_data` function.
#' 
#' @param msm.stan.mod The compiled Stan model to be used for each simulated 
#' data set. Should be an environment returned by `cmdstanr::cmdstan_model`
#' 
#' @param ... Additional arguments passed down to `cmdstanr::sample`
#'
#' @return Returns a tibble containing the summarized draws for each
#' of the models in the simulation
#' 
#' @export sim_msm_bayes_outcome
#'
sim_msm_bayes_outcome <- function(msm.data, msm.stan.mod, ...) {
  
  # Fit the Outcome-Stage Model
  msm_sim_fit <- msm_stan_mod$sample(
    data = stan.data,
    refresh = 0,
    sig_figs = 5,
    parallel_chains = 4,
    chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    max_treedepth = 12,
    adapt_delta = 0.9,
    show_messages = FALSE,
    ...
  )
  
  # Calculate a summary of the draws
  msm_result <- posterior::summarise_draws(
    msm_sim_fit$draws(
      variables = c("lp__", "b_Intercept", "b", "w_tilde", "Intercept", "sigma")
      )
    )
  
  # Return the update data frame
  return(msm_result)
}

