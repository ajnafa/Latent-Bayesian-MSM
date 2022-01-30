#-----------------------Project Data Helper Functions---------------------------
#-Author: A. Jordan Nafa-----------------------------Created: January 27, 2022-#
#-R Version: 4.1.0-----------------------------------Revised: January 27, 2022-#

# Custom theme for data visualizations
plot_theme <- function(...) {
  theme_bw() + theme(
    # Set the outer margins of the plot to 1/5 of an inch on all sides
    #plot.margin = margin(0.2, 0.2, 0.2, 0.2, "in"),
    # Specify the default settings for the plot title
    plot.title = element_text(
      size = 30,
      face = "bold",
      family = "serif"
    ),
    # Specify the default settings for caption text
    plot.caption = element_text(
      size = 12,
      family = "serif"
    ),
    # Specify the default settings for subtitle text
    plot.subtitle = element_text(
      size = 20,
      family = "serif"
    ),
    # Specify the default settings for axis titles
    axis.title = element_text(
      size = 22,
      face = "bold",
      family = "serif"
    ),
    # Specify the default settings specific to the x axis title
    axis.title.y = element_text(margin = margin(r = 10, l = -10)),
    # Specify the default settings specific to the y axis title
    axis.title.x = element_text(margin = margin(t = 10, b = -10)),
    # Specify the default settings for x axis text
    axis.text.x = element_text(
      size = 12,
      family = "serif"
    ),
    # Specify the default settings for y axis text
    axis.text.y = element_text(
      size = 12,
      family = "serif"
    ),
    # Specify the default settings for legend titles
    legend.title = element_text(
      size = 16,
      face = "bold",
      family = "serif"
    ),
    # Specify the default settings for legend text
    legend.text = element_text(
      size = 14,
      family = "serif"
    ),
    # Additional Settings Passed to theme()
    ...
  )
}

# A Function  for creating a field containing the meta data from a brms object----
stan_metadata <- function(x, ...){
  # Construct a field with relevant metadata
  meta_data <- map_dfr(
    .x = x$fit@stan_args,
    .f = ~ tibble(
      warmup = str_remove_all(.x$time_info[1], "[^[:digit:]|\\.]"),
      sampling = str_remove_all(.x$time_info[2], "[^[:digit:]|\\.]"),
      total = str_remove_all(.x$time_info[3], "[^[:digit:]|\\.]"),
      misc = str_remove_all(.x$time_info[4], "[^[:digit:]|\\.]"),
      metadata = c(
        str_c("stanc_version:", .x$stanc_version[1], sep = " "),
        str_c("opencl_device_name:", .x$opencl_device_name[1], sep = " "),
        str_c("opencl_platform_name", .x$opencl_platform_name[1], sep = " "),
        str_c("date", .x$start_datetime[1], sep = " ")
      )
    ),
    .id = "chain"
  )
  # Return the metadata field
  return(meta_data)
}

# A function for recovering the prior scales based on rstanarm's approach
prior_scaling <- function(y, df, ...) {
  # Initialize a matrix to store the location and scale values in
  scales <- matrix(NA, nrow = dim(df)[2] + 1, ncol = 2)
  y_mean <- mean(y, na.rm = TRUE) # Recover the mean of the response
  y_sd <- sd(y, na.rm = TRUE) # Recover the sd of the response
  
  # Priors for the Population-Level Intercept
  scales[1, ] <- c(y_mean, y_sd*2)
  
  # Priors for the random effects SDs
  scales[(dim(df)[2] + 1), 2] <- 1/y_sd
  
  # Retrieve the prior scales for the coefficients
  for (i in 2:dim(df)[2]) {
    scales[i, ] <- c(0.00, if_else(
      diff(range(df[, i])) > 1,
      (y_sd/sd(df[, i], na.rm = TRUE))*2,
      (1/sd(df[, i]))*2,
    ))
  }
  # Return the matrix object
  return(scales)
}
