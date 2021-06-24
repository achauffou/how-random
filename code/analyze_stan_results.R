# Generic function to analyse Stan results =====================================
#' Determine and use the correct function to analyse a Stan outcome
#'
analyse_stan_res <- function(spec, data, start, fits) {
  # Determine results folder and read RDS files:
  res_folder <- dirname(fits[[1]])
  cmdstan_fit <- readRDS(fits[[1]])
  rstan_fit <- readRDS(fits[[2]])
  data <- readRDS(data)
  start <-  readRDS(start)
  
  # Analyse simulation outcome:
  if (!file.exists(file.path(res_folder, "last_analysis.txt"))) {
    fun_name <- paste0("analyse_stan_res.", spec$stan_model)
    if (exists(fun_name)) {
      fun_name %>%
        get() %>%
        do.call(args = list(spec, data, start, cmdstan_fit, rstan_fit, res_folder))
      con <- file(file.path(res_folder, "last_analysis.txt"))
      writeLines(as.character(Sys.time()), con)
      close(con)
    }
  } else {
    message(paste("Results of model", spec$name, 
                  "seem to already have been analysed.",
                  "Skipping Stan results analysis..."))
  }
  
  # Return the current system time:
  as.character(Sys.time())
}


# Miscellaneous functions for the Stan simulations analyses ====================
#' Plot and save the posterior distribution and true value of several parameters
#'
stan_analyses_plot_save_params_post <- function(fit, params, res_folder) {
  lapply(params, function(x) {suppressWarnings({suppressMessages({
    suppressMessages({rstan::plot(fit, pars = x)}) %>%
      ggsave(file.path(res_folder, paste0("param_post_", x, ".pdf")), .,
             device = "pdf")
  })})})
}


# Functions to analyse specific Stan simulations ===============================
#' Analyse pollination binomial with intercepts only
#'
analyse_stan_res.pol_binom_02 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "beta", "gamma_pla", "gamma_pol", "sigma_beta", "sigma_gamma_pla", 
    "sigma_gamma_pol") %>%
    stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
  
  # Compute link:
  link <- compute_save_link(data, rstan_fit, link.pol_binom_02, res_folder)
  
  # Compute and plot AUC/ROC:
  stan_analyses_auc(data$Y_array$Y, link, res_folder, nb_samples = 100)
  stan_analyses_roc(data$Y_array$Y, link, res_folder, nb_samples = 100)
}

#' Analyse pollination binomial with intercepts and slopes
#'
analyse_stan_res.pol_binom_03 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda_bar", "beta", "gamma_pla", "gamma_pol",
    "lambda", "sigma_beta", "sigma_gamma_pla", "sigma_gamma_pol",
    "sigma_lambda") %>%
    stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
  
  # Compute link:
  link <- compute_save_link(data, rstan_fit, link.pol_binom_03, res_folder)
  
  # Compute and plot AUC/ROC:
  stan_analyses_auc(data$Y_array$Y, link, res_folder, nb_samples = 100)
  stan_analyses_roc(data$Y_array$Y, link, res_folder, nb_samples = 100)
}

#' Analyse pollination binomial with single lambda for all sites
#'
analyse_stan_res.pol_binom_04 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda", "beta", "gamma_pla", "gamma_pol",
    "sigma_beta", "sigma_gamma_pla", "sigma_gamma_pol") %>%
    stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
  
  # Compute link:
  link <- compute_save_link(data, rstan_fit, link.pol_binom_04, res_folder)
  
  # Compute and plot AUC/ROC:
  stan_analyses_auc(data$Y_array$Y, link, res_folder, nb_samples = 100)
  stan_analyses_roc(data$Y_array$Y, link, res_folder, nb_samples = 100)
}

#' Analyse pollination binomial with alpha, betas and lambdas
#'
analyse_stan_res.pol_binom_05 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda_bar", "beta", "lambda", "sigma_beta", "sigma_lambda") %>%
    stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
  
  # Compute link:
  link <- compute_save_link(data, rstan_fit, link.pol_binom_05, res_folder)
  
  # Compute and plot AUC/ROC:
  stan_analyses_auc(data$Y_array$Y, link, res_folder, nb_samples = 100)
  stan_analyses_roc(data$Y_array$Y, link, res_folder, nb_samples = 100)
}

#' Analyse pollination binomial with alpha, betas and single lambda
#'
analyse_stan_res.pol_binom_06 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda", "beta", "sigma_beta") %>%
    stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
  
  # Compute link:
  link <- compute_save_link(data, rstan_fit, link.pol_binom_06, res_folder)
  
  # Compute and plot AUC/ROC:
  stan_analyses_auc(data$Y_array$Y, link, res_folder, nb_samples = 100)
  stan_analyses_roc(data$Y_array$Y, link, res_folder, nb_samples = 100)
}
