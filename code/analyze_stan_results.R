# Generic function to analyse Stan results =====================================
#' Determine and use the correct function to analyse a Stan outcome
#'
analyse_stan_res <- function(spec, data, start, fits) {
  # Determine results folder:
  res_folder <- dirname(fits[[1]])
  
  # Analyse simulation outcome:
  if (!file.exists(file.path(res_folder, "last_analysis.txt"))) {
    fun_name <- paste0("analyse_stan_res.", spec$stan_model)
    if (exists(fun_name)) {
      # Read RDS files:
      cmdstan_fit <- readRDS(fits[[1]])
      rstan_fit <- readRDS(fits[[2]])
      data <- readRDS(data)
      start <-  readRDS(start)
      
      # Analyse results with the appropriate method:
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

#' Compute the variance of the predicted values of a Bayesian Bernoulli model
#' 
bayes_var_fit <- function(ypred) {
  apply(ypred, 1, var)
}

#' Compute the residual variance of a Bayesian Bernoulli model
#' 
bayes_var_res <- function(ypred) {
  (1 / ncol(ypred)) * apply(ypred, 1, function(x) sum(x * (1 - x)))
}

#' Compute the Bayesian R-squared from the modelled and residual variance
#' 
bayes_R2 <- function(var_fit, var_res) {
  var_fit / (var_fit + var_res)
}

#' Compute the Bayesian R2 statistics of a model (by groups)
#' 
calc_bayes_R2_stats <- function(fit, Y_array, groups, names) {
  # Get all link values:
  ypred <- rstan::extract(fit, pars = "link")[[1]]
  
  # Compute and save stats for the overall dataset:
  var_fit <- bayes_var_fit(ypred)
  var_res <- bayes_var_res(ypred)
  bayes_R2 <- bayes_R2(var_fit, var_res)
  R2 <- rbind(
    data.table::data.table(
      stat = "var_fit", group = "all", id = NA_integer_, name = NA_character_, 
      value = var_fit
    ),
    data.table::data.table(
      stat = "var_res", group = "all", id = NA_integer_, name = NA_character_, 
      value = var_res
    ),
    data.table::data.table(
      stat = "bayes_R2", group = "all", id = NA_integer_, name = NA_character_, 
      value = bayes_R2
    )
  )
  
  # Compute statistics group by group:
  for (the_group in groups) {
    the_names <- names[[the_group]]
    ids <- 1:length(the_names)
    nb_cores <- get_nb_cpus()
    new_R2 <- parallel::mclapply(ids, function(the_id) {
      the_name <- the_names[the_id]
      the_ypred <- ypred[, which(Y_array[[the_group]] == the_id)]
      the_var_fit <- bayes_var_fit(the_ypred)
      the_var_res <- bayes_var_res(the_ypred)
      the_bayes_R2 <- bayes_R2(the_var_fit, the_var_res)
      rbind(
        data.table::data.table(
          stat = "var_fit", group = the_group, id = the_id, name = the_name, 
          value = the_var_fit
        ),
        data.table::data.table(
          stat = "var_res", group = the_group, id = the_id, name = the_name, 
          value = the_var_res
        ),
        data.table::data.table(
          stat = "bayes_R2", group = the_group, id = the_id, name = the_name, 
          value = the_bayes_R2
        )
      )
    }, mc.cores = nb_cores) %>% data.table::rbindlist()
    R2 %<>% rbind(new_R2)
  }
  R2
}

#' Compute the WAIC of a model
#' 
calc_waic <- function(fit) {
  loo::extract_log_lik(fit, merge_chains = FALSE) %>%
    loo::waic()
}

#' Compute the LOO-IC of a model
#' 
calc_looic <- function(fit) {
  nb_cores <- get_nb_cpus()
  loglik <- loo::extract_log_lik(fit, merge_chains = FALSE)
  r_eff <- loo::relative_eff(exp(loglik), cores = nb_cores)
  loo::loo(loglik, r_eff = r_eff, cores = nb_cores)
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
  
  # Compute and save the R-squared statistics by group:
  calc_bayes_R2_stats(rstan_fit, data$Y_array, c("pla_id", "pol_id", "site_id"), list(
    pla_id = data$pla_names, pol_id = data$pol_names, site_id = data$site_names
  )) %>% saveRDS(file = file.path(res_folder, "bayes_R2_stats.rds"))
  
  # Compute and save WAIC:
  calc_waic(rstan_fit) %>% 
    saveRDS(file = file.path(res_folder, "waic.rds"))
  
  # Compute and save WAIC:
  calc_looic(rstan_fit) %>% 
    saveRDS(file = file.path(res_folder, "looic.rds"))
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
  
  # Compute and save the R-squared statistics by group:
  calc_bayes_R2_stats(rstan_fit, data$Y_array, c("pla_id", "pol_id", "site_id"), list(
    pla_id = data$pla_names, pol_id = data$pol_names, site_id = data$site_names
  )) %>% saveRDS(file = file.path(res_folder, "bayes_R2_stats.rds"))
  
  # Compute and save WAIC:
  calc_waic(rstan_fit) %>% 
    saveRDS(file = file.path(res_folder, "waic.rds"))
  
  # Compute and save WAIC:
  calc_looic(rstan_fit) %>% 
    saveRDS(file = file.path(res_folder, "looic.rds"))
}

#' Analyse pollination binomial with alpha, betas and lambdas
#'
analyse_stan_res.pol_binom_05 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda_bar", "beta", "lambda", "sigma_beta", "sigma_lambda") %>%
    stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
  
  # Compute and save the R-squared statistics by group:
  calc_bayes_R2_stats(rstan_fit, data$Y_array, c("pla_id", "pol_id", "site_id"), list(
    pla_id = data$pla_names, pol_id = data$pol_names, site_id = data$site_names
  )) %>% saveRDS(file = file.path(res_folder, "bayes_R2_stats.rds"))
  
  # Compute and save WAIC:
  calc_waic(rstan_fit) %>% 
    saveRDS(file = file.path(res_folder, "waic.rds"))
  
  # Compute and save WAIC:
  calc_looic(rstan_fit) %>% 
    saveRDS(file = file.path(res_folder, "looic.rds"))
}

#' Analyse pollination binomial with alpha, betas and single lambda
#'
analyse_stan_res.pol_binom_06 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda", "beta", "sigma_beta") %>%
    stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
  
  # Compute and save the R-squared statistics by group:
  calc_bayes_R2_stats(rstan_fit, data$Y_array, c("pla_id", "pol_id", "site_id"), list(
    pla_id = data$pla_names, pol_id = data$pol_names, site_id = data$site_names
  )) %>% saveRDS(file = file.path(res_folder, "bayes_R2_stats.rds"))
  
  # Compute and save WAIC:
  calc_waic(rstan_fit) %>% 
    saveRDS(file = file.path(res_folder, "waic.rds"))
  
  # Compute and save WAIC:
  calc_looic(rstan_fit) %>% 
    saveRDS(file = file.path(res_folder, "looic.rds"))
}


#' Analyse pollination binomial with alpha, betas, gammas product and lambdas
#'
analyse_stan_res.pol_binom_07 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda_bar", "beta", "zgamma_pla", "zgamma_pol",
    "lambda", "sigma_beta", "sigma_gamma", "sigma_lambda") %>%
    stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
  
  # Compute and save the R-squared statistics by group:
  calc_bayes_R2_stats(rstan_fit, data$Y_array, c("pla_id", "pol_id", "site_id"), list(
    pla_id = data$pla_names, pol_id = data$pol_names, site_id = data$site_names
  )) %>% saveRDS(file = file.path(res_folder, "bayes_R2_stats.rds"))
  
  # Compute and save WAIC:
  calc_waic(rstan_fit) %>% 
    saveRDS(file = file.path(res_folder, "waic.rds"))
  
  # Compute and save WAIC:
  calc_looic(rstan_fit) %>% 
    saveRDS(file = file.path(res_folder, "looic.rds"))
}

#' Analyse pollination binomial with alpha, betas, gammas prod and single lambda
#'
analyse_stan_res.pol_binom_08 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda", "beta", "zgamma_pla", "zgamma_pol",
    "sigma_beta", "sigma_gamma") %>%
    stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
  
  # Compute and save the R-squared statistics by group:
  calc_bayes_R2_stats(rstan_fit, data$Y_array, c("pla_id", "pol_id", "site_id"), list(
    pla_id = data$pla_names, pol_id = data$pol_names, site_id = data$site_names
  )) %>% saveRDS(file = file.path(res_folder, "bayes_R2_stats.rds"))
  
  # Compute and save WAIC:
  calc_waic(rstan_fit) %>% 
    saveRDS(file = file.path(res_folder, "waic.rds"))
  
  # Compute and save WAIC:
  calc_looic(rstan_fit) %>% 
    saveRDS(file = file.path(res_folder, "looic.rds"))
}
