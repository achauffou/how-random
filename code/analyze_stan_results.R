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
  compute_link <- TRUE
  rstan_file <- file.path(res_folder, "rstan-fit.rds")
  link_file <- file.path(res_folder, "link.rds")
  if (file.exists(link_file)) {
    if (file.info(link_file)$ctime > file.info(rstan_file)$ctime) {
      compute_link <- FALSE
      link <- readRDS(link_file)
    }
  }
  if (compute_link == TRUE) {
    link <- link.pol_binom_02(data, rstan_fit)
    saveRDS(link, link_file)
  }
  
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
  compute_link <- TRUE
  rstan_file <- file.path(res_folder, "rstan-fit.rds")
  link_file <- file.path(res_folder, "link.rds")
  if (file.exists(link_file)) {
    if (file.info(link_file)$ctime > file.info(rstan_file)$ctime) {
      compute_link <- FALSE
      link <- readRDS(link_file)
    }
  }
  if (compute_link == TRUE) {
    link <- link.pol_binom_03(data, rstan_fit)
    saveRDS(link, link_file)
  }
  
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
  compute_link <- TRUE
  rstan_file <- file.path(res_folder, "rstan-fit.rds")
  link_file <- file.path(res_folder, "link.rds")
  if (file.exists(link_file)) {
    if (file.info(link_file)$ctime > file.info(rstan_file)$ctime) {
      compute_link <- FALSE
      link <- readRDS(link_file)
    }
  }
  if (compute_link == TRUE) {
    link <- link.pol_binom_04(data, rstan_fit)
    saveRDS(link, link_file)
  }
  
  # Compute and plot AUC/ROC:
  stan_analyses_auc(data$Y_array$Y, link, res_folder, nb_samples = 100)
  stan_analyses_roc(data$Y_array$Y, link, res_folder, nb_samples = 100)
}

#' Link function of pollination binomial with single lambda for all sites
#'
link.pol_binom_04 <- function(data, fit) {
  alpha <- as.matrix(fit, pars = "alpha")
  beta <- as.matrix(fit, pars = "beta")
  gamma_pla <- as.matrix(fit, pars = "gamma_pla")
  gamma_pol <- as.matrix(fit, pars = "gamma_pol")
  lambda <- as.matrix(fit, pars = "lambda")
  parallel::mclapply(1:(dim(fit)[1] * dim(fit)[2]), function(x) {
    p <- boot::inv.logit(
      alpha[x] + beta[x, data$Y_array$site_id] +
        gamma_pla[x, data$Y_array$pla_id] + gamma_pol[x, data$Y_array$pol_id] +
        lambda * data$SS
    )
    names(p) <- NULL
    p
  }, mc.cores = get_nb_cpus()) %>% do.call(rbind, args = .)
}
