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
  lapply(params, function(x) {suppressMessages({
    stan_sim_analyses_plot_param_post_true(fit, data, x) %>%
      ggsave(file.path(res_folder, paste0("param_post_true_ ", x, ".pdf")), ., 
             device = "pdf")
  })})
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
}
