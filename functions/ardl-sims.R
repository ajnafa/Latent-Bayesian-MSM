#' Function for Fitting the ARDL Model Simulations
#' 
#' @param ardl.data A list object containing the data to be passed to the
#' Stan model as returned by the `make_ardl_data` function.
#' 
#' @param ardl.mod The compiled Stan model to be used for each simulated data
#' set. Should be an environment returned by `cmdstanr::cmdstan_model`
#' 
#' @param ... Additional arguments passed down to `cmdstanr::sample`
#'
#' @return Returns a tibble containing the summarized draws for each
#' of the models in the simulation
#' 
#' @export sim_ardl_bayes
#'
sim_ardl_bayes <- function(ardl.data, ardl.mod, ...) {
  
  # Fit the Outcome-Stage Model
  ardl_fit <- ardl.mod$sample(
    data = ardl.data,
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
  ardl_result <- summarise_draws(ardl_fit$draws())
  
  # Return the data frame of draws
  return(ardl_result)
}

#' Function for Building the Data to Pass to the ARDL Stan Model
#'
#' @param sim.data A simulated data set returned by the `dgp_sim` function
#' 
#' @param ardl.form A formula for the ARDL model to be passed down to 
#' `brms::make_standata`
#' 
#' @param prior.scale Scale factor to be used for auto-scaling of the priors
#' on the coefficients. For details see `vignette("priors", package = "rstanarm")`
#' 
#' @param ... Currently unused
#'
#' @return Returns a list object with the data to be passed to the Stan model
#' 
#' @export make_ardl_data
#'
make_ardl_data <- function(sim.data, ardl.form, prior.scale, ...) {
  # Take advantage of brms functionality to get the initial data
  ardl_data <- brms::make_standata(
    ardl.form,
    family = brms::brmsfamily(family = "gaussian", link = "identity"), 
    data = sim.data
  )
  
  # Priors on the coefficients
  beta_sds <- prior.scale * (sd(ardl_data$Y)/apply(ardl_data$X, 2, sd))
  
  # Prepare the data for use with Stan
  ardl_data <- list(
    N = ardl_data$N,
    K = ardl_data$K,
    Y = ardl_data$Y,
    X = ardl_data$X,
    beta_sds = beta_sds
  )
  
  # Return the list of data
  return(ardl_data)
}
