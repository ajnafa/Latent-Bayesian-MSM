#' Function for Simulating the Datasets, Modified Version of Blackwell and Glynn (2018)
#'
#' @param .groups Number of groups in the cross-sectional dimension of the data
#' 
#' @param .periods Number of periods in the temporal dimension of the data
#' 
#' @param .true_gamma Numeric argument(s) indicating true value of the lagged treatment
#' effect.
#' 
#' @param .treat_conf Logical argument(s) indicating whether the unobserved confounder 
#' is correlated with the treatment
#' 
#' @param .seed Interger value to use for the random number seed
#' 
#' @param ... Additional arguments for the simulation, currently unused
#'
#' @export dgp_sim
#'
dgp_sim <- function(.groups, .periods, .true_gamma, .treat_conf, .seed, ...) { 
  
  # Requires data.table package
  require(data.table)
  
  # Set the rng seed
  set.seed(.seed)
  
  # Define the fixed parameter values
  mu = c(0, -0.1, -0.1) # Effect of X[t-1], X[t], and X[t] * X[t-1] on Y[t]
  alpha_t0 = c(-0.23, 2.5) # Values of X at time t = 1
  alpha_t1 = c(-1.3, 2.5, 1.5) # Values of X at time t = 2:T
  gamma = .true_gamma
  
  # Create N x T Matricies
  Y = X = Z = matrix(NA_real_, nrow = .groups, ncol = .periods)
  
  # Country-Specific Time Invariant Confounding
  U = rnorm(.groups, sd = 0.1)
  upsilon = .9 * U
  
  # Simulate the first time period for each country
  Y_00 = 0.8 + upsilon + rnorm(.groups, sd = 0.1)
  Y_01 = 0.8 + mu[2] + upsilon + rnorm(.groups, sd = 0.1)
  Y_10 = Y_11 = Y_00
  
  # Z[t[1]] ~ 0.8 + upsilon + N(0.4, 0.1)
  Z[, 1] = 0.8 + upsilon + rnorm(.groups, 0.4, 0.1)
  
  # Prob of X[t[1]] = 1 with time-invariant confounding
  if (isTRUE(.treat_conf)) {
    X_prob = plogis(alpha_t0[1] + upsilon + alpha_t0[2] * Z)
    X[, 1] = rbinom(.groups, 1, prob = X_prob)
  } 
  
  # Prob of X[t[1]] = 1 without time-invariant confounding
  else {
    X_prob = plogis(alpha_t0[1] + alpha_t0[2] * Z)
    X[, 1] = rbinom(.groups, 1, prob = X_prob)
  }
  
  # Actual Values of Y[t[1]]
  Y[, 1] = Y_01 * X[, 1] + Y_00 * (1 - X[, 1])
  
  # Recursively simulate the data for time t = 2:T
  for (i in 2:.periods) {
    
    # Potential Outcomes Y(a, a')
    Y_11 = 0.8 + upsilon + mu[1] + mu[3] + rnorm(.groups, sd = 0.1)
    Y_10 = 0.8 + upsilon + mu[1] + rnorm(.groups, sd = 0.1)
    Y_01 = 0.8 + upsilon + mu[2] + rnorm(.groups, sd = 0.1)
    Y_00 = 0.8 + upsilon + rnorm(.groups, sd = 0.1)
    
    # Potential Confounders Z(a, a')
    Z_1 = 1 * gamma + 0.5 + 0.7 * U + rnorm(.groups, sd = 0.1)
    Z_0 = 0 * gamma + 0.5 + 0.7 * U + rnorm(.groups, sd = 0.1)
    Z[, i] = X[, i - 1] * Z_1 + (1 - X[, i - 1]) * Z_0
    
    # Treatment X with time-invariant confounding
    if (isTRUE(.treat_conf)) {
      X_pr = alpha_t1[1] + alpha_t1[2] * Z[, i] + alpha_t1[3] * Y[, i - 1] + upsilon + rnorm(.groups, 0, 1)
      X[, i] = 1 * (X_pr > 0)
    }
    
    # Treatment X without time-invariant confounding
    else {
      X_pr = alpha_t1[1] + alpha_t1[2] * Z[, i] + alpha_t1[3] * Y[, i - 1] + rnorm(.groups, 0, 1)
      X[, i] = 1 * (X_pr > 0)
    }
    
    # The Response Vector Y
    Y[, i] = # Control
      Y_00 +
      # Effect of Lagged Treatment X
      X[, i - 1] * (Y_10 - Y_00) +
      # Effect of Contemporaneous X
      X[, i] * (Y_01 - Y_00) + 
      # Effect of Contemporaneous X and Lagged X
      X[, i - 1] * X[, i] * ((Y_11 - Y_01) - (Y_10 - Y_00)) 
  }
  
  # Reshape and Combine the Matrices
  out = list(X, Y, Z)
  
  for (i in 1:3) {
    out[[i]] = as.data.table(out[[i]])
    out[[i]][, unit := 1:.N]
    out[[i]] = melt(out[[i]], id.vars = "unit", variable.name = "time")
    out[[i]] = out[[i]][, time := as.numeric(gsub("V", "", time))]
  }
  
  colnames(out[[1]])[3] = "X"
  colnames(out[[2]])[3] = "Y"
  colnames(out[[3]])[3] = "Z"
  
  out = Reduce(function(x, y) merge(x, y, by = c("unit", "time")), out)
  
  # Generate the Lags
  out[, `:=` (
    Y_Lag = shift(Y),
    Y_Lag_2 = shift(Y, 2),
    X_Lag = shift(X),
    Z_Lag = shift(Z))]
  
  # Exclude missing data due to lags
  out = out[!is.na(Y_Lag_2)]
  
  return(out)
}

