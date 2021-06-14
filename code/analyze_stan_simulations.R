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
    ylim(-2.5, 2.5) +
    theme_classic()
  if (length(unique(nb_data)) == 1) {
    plt <- plt + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                axis.ticks.x = element_blank(), axis.line.x = element_blank())
  }
  plt
}

#' Plot the posterior distribution and true value of a parameter
#' 
misc_stan_analyses_plot_param_post_true <- function(fit, data, param) {
  true_vals <- data[[param]]
  suppressWarnings({rstan::plot(fit, pars = param) +
    geom_point(
      aes(x=x, y=y), data.frame(y = length(true_vals):1, x = true_vals), 
      color = "blue") +
    xlim(-2.5, 2.5)})
}

#' Plot and save the posterior distribution and true value of several parameters
#' 
misc_stan_analyses_plot_save_params_post_true <- function(
  fit, data, params, res_folder
) {
  suppressWarnings({lapply(params, function(x) {
    misc_stan_analyses_plot_param_post_true(fit, data, x) %>%
      ggsave(file.path(res_folder, paste0("param_post_true_ ", x, ".pdf")), ., 
             device = "pdf")
  })})
}


# Functions to analyse specific Stan simulations ===============================
#' Analyse simulation results of Stan model pol_logit_f
#' 
analyse_stan_sim.pol_logit_f <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Perform same analyses as pol_logit_h:
  analyse_stan_sim.pol_logit_h(
    spec, data, start, cmdstan_fit, rstan_fit, res_folder
  )
  
  # Plot the error distribution of parameters lambda and nu:
  misc_stan_analyses_plot_param_error(
    rstan_fit, "lambda", data$lambda, 
    data$Y_array[, .N, by = .(site_id)][order(site_id)][['N']], FALSE) %>%
    ggsave(file.path(res_folder, "param_errors_lambda.pdf"), ., device = "pdf")
  misc_stan_analyses_plot_param_error(
    rstan_fit, "nu", data$nu, 
    data$Y_array[, .N, by = .(site_id)][order(site_id)][['N']], FALSE) %>%
    ggsave(file.path(res_folder, "param_errors_nu.pdf"), ., device = "pdf")
  
  # Plot the posterior distribution and true value of lambda and nu:
  misc_stan_analyses_plot_save_params_post_true(
    rstan_fit, data, 
    c("lambda", "nu", "sigma_lambda", "sigma_nu"), 
    res_folder
  )
}

#' Analyse simulation results of Stan model pol_logit_g
#'
analyse_stan_sim.pol_logit_g <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Perform same analyses as pol_logit_f:
  analyse_stan_sim.pol_logit_f(
    spec, data, start, cmdstan_fit, rstan_fit, res_folder
  )
  
  # Plot the error distribution of the parameter alpha:
  misc_stan_analyses_plot_param_error(
    rstan_fit, "alpha", data$alpha, 
    nrow(data$Y_array), FALSE) %>%
    ggsave(file.path(res_folder, "param_errors_alpha.pdf"), ., device = "pdf")
  
  # Plot the posterior distribution and true value of alpha:
  misc_stan_analyses_plot_save_params_post_true(
    rstan_fit, data, "alpha", res_folder
  )
}

#' Analyse simulation results of Stan model pol_logit_h
#'
analyse_stan_sim.pol_logit_h <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot error distributions of beta, gamma_pla, gamma_pol
  misc_stan_analyses_plot_param_error(
    rstan_fit, "beta", data$alpha, 
    nrow(data$Y_array), FALSE) %>%
    ggsave(file.path(res_folder, "param_errors_alpha.pdf"), ., device = "pdf")
  misc_stan_analyses_plot_param_error(
    rstan_fit, "beta", data$beta, 
    data$Y_array[, .N, by = .(site_id)][order(site_id)][['N']], FALSE) %>%
    ggsave(file.path(res_folder, "param_errors_beta.pdf"), ., device = "pdf")
  misc_stan_analyses_plot_param_error(
    rstan_fit, "gamma_pla", data$gamma_pla, 
    data$Y_array[, .N, by = .(pla_id)][order(pla_id)][['N']], FALSE) %>%
    ggsave(file.path(res_folder, "param_errors_gamma_pla.pdf"), ., device = "pdf")
  misc_stan_analyses_plot_param_error(
    rstan_fit, "gamma_pol", data$gamma_pol, 
    data$Y_array[, .N, by = .(pol_id)][order(pol_id)][['N']], FALSE) %>%
    ggsave(file.path(res_folder, "param_errors_gamma_pol.pdf"), ., device = "pdf")
  
  # Plot the posterior distribution and true value of lambda and nu:
  misc_stan_analyses_plot_save_params_post_true(
    rstan_fit, data, 
    c("alpha", "beta", "sigma_beta", "gamma_pla", "sigma_gamma_pla", 
      "gamma_pol", "sigma_gamma_pol"), 
    res_folder
  )
}

#' Analyse simulation results of Stan model pol_logit_i
#'
analyse_stan_sim.pol_logit_i <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Perform same analyses as pol_logit_g:
  analyse_stan_sim.pol_logit_g(
    spec, data, start, cmdstan_fit, rstan_fit, res_folder
  )
  
  # Plot the posterior distribution and true value of lambda_bar and nu_bar:
  misc_stan_analyses_plot_save_params_post_true(
    rstan_fit, data, 
    c("lambda_bar", "nu_bar"), 
    res_folder
  )
}
