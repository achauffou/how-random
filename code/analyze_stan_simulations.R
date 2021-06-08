# Generic function to analyse Stan simulations =================================
#' determine and use the correct function to analyse a Stan simulation outcome
#' 
analyse_stan_sim <- function(spec, data, start, fits) {
  # Determine results folder and read RDS files:
  res_folder <- dirname(fits[[1]])
  cmdstan_fit <- readRDS(fits[[1]])
  rstan_fit <- readRDS(fits[[2]])
  data <- readRDS(data)
  start <-  readRDS(start)
  
  # Analyse simulation outcome:
  fun_name <- paste0("analyse_stan_sim.", spec$stan_model)
  if (exists(fun_name)) {
    fun_name %>%
      get() %>%
      do.call(args = list(spec, data, start, cmdstan_fit, rstan_fit, res_folder))
  }
  
  # Return RStan fit summary:
  list(rstan::summary(rstan_fit)$summary)
}


# Miscellaneous functions for the Stan simulations analyses ====================
#' Plot the errors of a parameter ordered by its number of data (mean + 95% CI)
#' 
misc_stan_analyses_plot_param_error <- function(
  fit, param, exact_values, nb_data, aggregate_by_nb_data = TRUE
) {
  # Create the data.table for the plot:
  plt_data <- as.data.frame(fit, pars = param) %>%
    data.table::as.data.table() %>%
    data.table::melt(id.vars = character())
  
  # Add variable ID, error and number of data:
  if (length(exact_values) > 1) {
    plt_data[, id := as.integer(stringr::str_extract(variable, "[0-9]+"))]
  } else {
    plt_data[, id := 1]
  }
  plt_data[, ':='(
    error = value - exact_values[id],
    nb_data = nb_data[id]
  )]
  
  # Make the boxplot:
  plt <- ggplot(plt_data, aes_string(
    x = "nb_data", y = "error", 
    group = ifelse(aggregate_by_nb_data, "nb_data", "id"), width = 0.4
  )) +
    stat_summary(geom = "point", fun = mean, position=position_dodge(.9)) +
    stat_summary(geom = "errorbar", fun.min = function(x) quantile(x, 0.05),
                 fun.max = function(x) quantile(x, 0.95), position=position_dodge(.9)) +
    xlab("Amount of data") +
    ylab("Error") +
    theme_classic()
  if (length(unique(nb_data)) == 1) {
    plt <- plt + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                axis.ticks.x = element_blank(), axis.line.x = element_blank())
  }
  plt
}
