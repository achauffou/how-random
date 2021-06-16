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
#' Plot the posterior distribution and true value of a parameter
#' 
stan_sim_analyses_plot_param_post_true <- function(
  fit, data, param
) {
  true_vals <- data[[param]]
  suppressMessages({rstan::plot(fit, pars = param) +
    geom_point(
      aes(x=x, y=y), data.frame(y = length(true_vals):1, x = true_vals), 
      color = "blue") +
    xlim(-3.0, 3.0)})
}

#' Plot and save the posterior distribution and true value of several parameters
#' 
stan_sim_analyses_plot_save_params_post_true <- function(
  fit, data, params, res_folder
) {
  lapply(params, function(x) {suppressWarnings({suppressMessages({
    stan_sim_analyses_plot_param_post_true(fit, data, x) %>%
      ggsave(file.path(res_folder, paste0("param_post_true_ ", x, ".pdf")), ., 
             device = "pdf")
  })})})
}

#' Plot the errors of a parameter ordered by its number of data (mean + 95% CI)
#' 
stan_sim_analyses_plot_param_error <- function(
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
    geom_hline(yintercept = 0.0, color = "grey") +
    stat_summary(geom = "point", fun = mean, position=position_dodge(.9)) +
    stat_summary(geom = "errorbar", fun.min = function(x) quantile(x, 0.05),
                 fun.max = function(x) quantile(x, 0.95), position=position_dodge(.9)) +
    xlab("Amount of data") +
    ylab("Error") +
    ylim(-3.0, 3.0) +
    theme_classic()
  if (length(unique(nb_data)) == 1) {
    plt <- plt + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(), axis.line.x = element_blank())
  }
  plt
}

#' Plot the errors of several parameters ordered by their number of data
#' 
stan_sim_analyses_plot_save_params_errors <- function(
  fit, data, params, id_names, res_folder
) {
  lapply(1:length(params), function(x) {suppressWarnings({suppressMessages({
    stan_sim_analyses_plot_param_error(
      fit, params[x], data[[params[x]]], 
      data$Y_array[, .N, by = c(id_names[x])][order(get(id_names[x]))][['N']]) %>%
      ggsave(file.path(res_folder, paste0("param_error_ ", params[x], ".pdf")), ., 
             device = "pdf")
  })})})
}


# Functions to analyse specific Stan simulations ===============================
#' Analyse pollination binomial with centered intercepts only
#' 
analyse_stan_sim.pol_binom_01 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("beta", "gamma_pla", "gamma_pol", "sigma_beta", "sigma_gamma_pla",
    "sigma_gamma_pol") %>%
    stan_sim_analyses_plot_save_params_post_true(rstan_fit, data, ., res_folder)
  
  # Plot posterior error of multilevel parameters:
  stan_sim_analyses_plot_save_params_errors(
    rstan_fit, data, c("beta", "gamma_pla", "gamma_pol"), 
    c("site_id", "pla_id", "pol_id"), res_folder
  )
}

#' Analyse pollination binomial with intercepts only
#' 
analyse_stan_sim.pol_binom_02 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "beta", "gamma_pla", "gamma_pol", "sigma_beta", "sigma_gamma_pla",
    "sigma_gamma_pol") %>%
    stan_sim_analyses_plot_save_params_post_true(rstan_fit, data, ., res_folder)
  
  # Plot posterior error of multilevel parameters:
  stan_sim_analyses_plot_save_params_errors(
    rstan_fit, data, c("beta", "gamma_pla", "gamma_pol"), 
    c("site_id", "pla_id", "pol_id"), res_folder
  )
}

#' Analyse pollination binomial with intercepts and slopes
#' 
analyse_stan_sim.pol_binom_03 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda_bar", "nu_bar", "beta", "gamma_pla", "gamma_pol", 
    "lambda", "nu", "sigma_beta", "sigma_gamma_pla", "sigma_gamma_pol",
    "sigma_lambda", "sigma_nu") %>%
    stan_sim_analyses_plot_save_params_post_true(rstan_fit, data, ., res_folder)
  
  # Plot posterior error of multilevel parameters:
  stan_sim_analyses_plot_save_params_errors(
    rstan_fit, data, c("beta", "gamma_pla", "gamma_pol", "lambda", "nu"), 
    c("site_id", "pla_id", "pol_id", "site_id", "site_id"), res_folder
  )
}
