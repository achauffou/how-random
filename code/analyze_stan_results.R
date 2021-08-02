# Generic function to analyse Stan results =====================================
#' Determine and use the correct function to analyse a Stan outcome
#'
analyse_stan_res <- function(spec, data, start, fits) {
  # Determine results folder:
  res_folder <- dirname(fits[[1]])
  
  # Check whether simulations can be skipped:
  skip_analyses <- FALSE
  if (file.exists(file.path(res_folder, "prev_analyses.txt"))) {
    prev_modules <- readr::read_lines(file.path(res_folder, "prev_analyses.txt"))
  } else {
    prev_modules <- NA_character_
  }
  obj_name <- paste0("stan_res_mods.", spec$stan_model)
  if (exists(obj_name)) {
    if (all(get(obj_name) %in% prev_modules)) {
      skip_analyses <- TRUE
      message(paste(
        "Results of model", spec$name, "seem to already have been analysed.",
        "Skipping Stan results analysis..."
      ))
    }
  }
  
  # Analyse simulation outcome:
  if (skip_analyses == FALSE) {
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
        do.call(args = list(spec, data, start, cmdstan_fit, rstan_fit, 
                            res_folder, prev_modules))
    }
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
      if (length(which(Y_array[[the_group]] == the_id)) > 1) {
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
      } else {
        data.table::data.table(
          stat = character(), group = character(), id = integer(), 
          name = character(), value = numeric()
        )
      }
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


# Functions to analyse pollination binomial Stan results =======================
#' Miscellaneous analyses for pollination binomial Stan results
#' 
analyse_stan_res.misc_pol_binom <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules, include_origin = 0
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")

  # Compute and save the R-squared statistics by group:
  if (!"bayes_R2" %in% prev_modules) {
    if (include_origin == 0) {
      calc_bayes_R2_stats(
        rstan_fit, data$Y_array, c("pla_id", "pol_id", "site_id"),
        list(pla_id = data$pla_names, pol_id = data$pol_names,
             site_id = data$site_names)
      ) %>% saveRDS(file = file.path(res_folder, "bayes_R2_stats.rds"))
    } else if (include_origin == 1) {
      the_array <- data$Y_array
      the_array[, ':='(
        sp1_invasive = sp1_invasive + 1
      )]
      calc_bayes_R2_stats(
        rstan_fit, data$Y_array, c("pla_id", "pol_id", "site_id", "sp1_invasive"),
        list(pla_id = data$pla_names, pol_id = data$pol_names,
             site_id = data$site_names, sp1_invasive = c("Native", "Introduced"))
      ) %>% saveRDS(file = file.path(res_folder, "bayes_R2_stats.rds"))
    } else {
      the_array <- data$Y_array
      the_array[, ':='(
        sp1_invasive = sp1_invasive + 1,
        sp2_invasive = sp2_invasive + 1,
        invasive = sp1_invasive + sp2_invasive + 1
      )]
      calc_bayes_R2_stats(
        rstan_fit, data$Y_array, c("pla_id", "pol_id", "site_id", "sp1_invasive",
                                   "sp2_invasive", "invasive"),
        list(pla_id = data$pla_names, pol_id = data$pol_names,
             site_id = data$site_names, sp1_invasive = c("Native", "Introduced"),
             sp2_invasive = c("Native", "Introduced"),
             invasive = c("None", "Only one", "Both"))
      ) %>% saveRDS(file = file.path(res_folder, "bayes_R2_stats.rds"))
    }
    readr::write_lines("bayes_R2", prev_path, append = TRUE)
  }

  # Compute and save WAIC:
  if (!"waic" %in% prev_modules) {
    calc_waic(rstan_fit) %>%
      saveRDS(file = file.path(res_folder, "waic.rds"))
    readr::write_lines("waic", prev_path, append = TRUE)
  }

  # Compute and save WAIC:
  if (!"looic" %in% prev_modules) {
    calc_looic(rstan_fit) %>%
      saveRDS(file = file.path(res_folder, "looic.rds"))
    readr::write_lines("looic", prev_path, append = TRUE)
  }
}
stan_res_mods.misc_pol_binom <- c("bayes_R2", "waic", "looic")

#' Analyse pollination binomial with intercepts only
#'
analyse_stan_res.pol_binom_02 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "beta", "gamma_pla", "gamma_pol", "sigma_beta", 
      "sigma_gamma_pla", "sigma_gamma_pol") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
}
stan_res_mods.pol_binom_02 <- c("post_param_plots")

#' Analyse pollination binomial with intercepts and slopes
#'
analyse_stan_res.pol_binom_03 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda_bar", "beta", "gamma_pla", "gamma_pol",
      "lambda", "sigma_beta", "sigma_gamma_pla", "sigma_gamma_pol",
      "sigma_lambda") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous pollination results analyses:
  analyse_stan_res.misc_pol_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.pol_binom_03 <- c("post_param_plots", stan_res_mods.misc_pol_binom)

#' Analyse pollination binomial with single lambda for all sites
#'
analyse_stan_res.pol_binom_04 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda", "beta", "gamma_pla", "gamma_pol",
      "sigma_beta", "sigma_gamma_pla", "sigma_gamma_pol") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous pollination results analyses:
  analyse_stan_res.misc_pol_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.pol_binom_04 <- c("post_param_plots", stan_res_mods.misc_pol_binom)

#' Analyse pollination binomial with alpha, betas and lambdas
#'
analyse_stan_res.pol_binom_05 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda_bar", "beta", "lambda", "sigma_beta", "sigma_lambda") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous pollination results analyses:
  analyse_stan_res.misc_pol_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.pol_binom_05 <- c("post_param_plots", stan_res_mods.misc_pol_binom)

#' Analyse pollination binomial with alpha, betas and single lambda
#'
analyse_stan_res.pol_binom_06 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda", "beta", "sigma_beta") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous pollination results analyses:
  analyse_stan_res.misc_pol_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.pol_binom_06 <- c("post_param_plots", stan_res_mods.misc_pol_binom)

#' Analyse pollination binomial with alpha, betas, gammas product and lambdas
#'
analyse_stan_res.pol_binom_07 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda_bar", "beta", "zgamma_pla", "zgamma_pol",
      "lambda", "sigma_beta", "sigma_gamma", "sigma_lambda") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous pollination results analyses:
  analyse_stan_res.misc_pol_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.pol_binom_07 <- c("post_param_plots", stan_res_mods.misc_pol_binom)

#' Analyse pollination binomial with alpha, betas, gammas prod and single lambda
#'
analyse_stan_res.pol_binom_08 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda", "beta", "zgamma_pla", "zgamma_pol",
      "sigma_beta", "sigma_gamma") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous pollination results analyses:
  analyse_stan_res.misc_pol_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.pol_binom_08 <- c("post_param_plots", stan_res_mods.misc_pol_binom)

#' Analyse pollination binomial with two terms for bioclimatic suitability
#'
analyse_stan_res.pol_binom_09 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda_pla", "lambda_pol", "beta", "gamma_pla", "gamma_pol",
      "sigma_beta", "sigma_gamma_pla", "sigma_gamma_pol") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous pollination results analyses:
  analyse_stan_res.misc_pol_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.pol_binom_09 <- c("post_param_plots", stan_res_mods.misc_pol_binom)

#' Analyse pollination binomial with single lambda for all sites (2 origins)
#'
analyse_stan_res.pol_binom_14 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda", "mu_pla", "beta", "gamma_pla", "gamma_pol",
      "sigma_beta", "sigma_gamma_pla", "sigma_gamma_pol") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous pollination results analyses:
  analyse_stan_res.misc_pol_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules, include_origin = 1)
}
stan_res_mods.pol_binom_14 <- c("post_param_plots", stan_res_mods.misc_pol_binom)

#' Analyse pollination binomial with two terms for bioclimatic suitability (1 origin)
#'
analyse_stan_res.pol_binom_19 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda_pla", "lambda_pol", "mu_pla", "beta", 
      "gamma_pla", "gamma_pol", "sigma_beta", "sigma_gamma_pla", 
      "sigma_gamma_pol") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous pollination results analyses:
  analyse_stan_res.misc_pol_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules, include_origin = 1)
}
stan_res_mods.pol_binom_19 <- c("post_param_plots", stan_res_mods.misc_pol_binom)

#' Analyse pollination binomial with single lambda for all sites (2 origins)
#'
analyse_stan_res.pol_binom_24 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda", "mu_pla", "mu_pol", "beta", "gamma_pla", "gamma_pol",
      "sigma_beta", "sigma_gamma_pla", "sigma_gamma_pol") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous pollination results analyses:
  analyse_stan_res.misc_pol_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules, include_origin = 2)
}
stan_res_mods.pol_binom_24 <- c("post_param_plots", stan_res_mods.misc_pol_binom)

#' Analyse pollination binomial with two terms for bioclimatic suitability (2 origins)
#'
analyse_stan_res.pol_binom_29 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda_pla", "lambda_pol", "mu_pla", "mu_pol", "beta", 
      "gamma_pla", "gamma_pol", "sigma_beta", "sigma_gamma_pla", 
      "sigma_gamma_pol") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous pollination results analyses:
  analyse_stan_res.misc_pol_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules, include_origin = 2)
}
stan_res_mods.pol_binom_29 <- c("post_param_plots", stan_res_mods.misc_pol_binom)


# Functions to analyse all interactions binomial Stan results ==================
#' Miscellaneous analyses for all interactions binomial Stan results
#' 
analyse_stan_res.misc_all_binom <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules, include_origin = 0
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Compute and save the R-squared statistics by group:
  if (!"bayes_R2" %in% prev_modules) {
    if (include_origin == 0) {
      the_array <- data$Y_array
      the_array[, ':='(
        site_type = data$site_type[site_id]
      )]
      calc_bayes_R2_stats(
        rstan_fit, the_array, c("sp1_id", "sp2_id", "site_id", "site_type"),
        list(sp1_id = data$sp_names[data$sp_group %% 2 == 1],
             sp2_id = data$sp_names[data$sp_group %% 2 == 0],
             site_id = data$site_names,
             site_type = data$type_names)
      ) %>% saveRDS(file = file.path(res_folder, "bayes_R2_stats.rds"))
    } else if (include_origin == 1) {
      the_array <- data$Y_array
      the_array[, ':='(
        site_type = data$site_type[site_id],
        sp1_invasive = sp1_invasive + 1
      )]
      calc_bayes_R2_stats(
        rstan_fit, the_array, c("sp1_id", "sp2_id", "site_id", "site_type", "sp1_invasive"),
        list(sp1_id = data$sp_names[data$sp_group %% 2 == 1],
             sp2_id = data$sp_names[data$sp_group %% 2 == 0],
             site_id = data$site_names,
             site_type = data$type_names,
             sp1_invasive = c("Native", "Introduced"))
      ) %>% saveRDS(file = file.path(res_folder, "bayes_R2_stats.rds"))
    } else {
      the_array <- data$Y_array
      the_array[, ':='(
        site_type = data$site_type[site_id],
        sp1_invasive = sp1_invasive + 1,
        sp2_invasive = sp2_invasive + 1,
        invasive = sp1_invasive + sp2_invasive + 1
      )]
      calc_bayes_R2_stats(
        rstan_fit, the_array, c("sp1_id", "sp2_id", "site_id", "site_type", 
                                "sp1_invasive", "sp2_invasive", "invasive"),
        list(sp1_id = data$sp_names[data$sp_group %% 2 == 1],
             sp2_id = data$sp_names[data$sp_group %% 2 == 0],
             site_id = data$site_names,
             site_type = data$type_names,
             sp1_invasive = c("Native", "Introduced"),
             sp2_invasive = c("Native", "Introduced"),
             invasive = c("None", "Only one", "Both"))
      ) %>% saveRDS(file = file.path(res_folder, "bayes_R2_stats.rds"))
    }
    readr::write_lines("bayes_R2", prev_path, append = TRUE)
  }
  
  # Compute and save WAIC:
  if (!"waic" %in% prev_modules) {
    calc_waic(rstan_fit) %>%
      saveRDS(file = file.path(res_folder, "waic.rds"))
    readr::write_lines("waic", prev_path, append = TRUE)
  }
  
  # Compute and save WAIC:
  if (!"looic" %in% prev_modules) {
    calc_looic(rstan_fit) %>%
      saveRDS(file = file.path(res_folder, "looic.rds"))
    readr::write_lines("looic", prev_path, append = TRUE)
  }
}
stan_res_mods.misc_all_binom <- c("bayes_R2", "waic", "looic")

#' Analyse all interactions binomial with intercepts and slopes
#'
analyse_stan_res.all_binom_03 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda_bar", "beta", "gamma", "lambda", "sigma_beta", 
      "sigma_gamma", "sigma_lambda") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous all interactions results analyses:
  analyse_stan_res.misc_all_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.all_binom_03 <- c("post_param_plots", stan_res_mods.misc_all_binom)

#' Analyse all interactions binomial with single lambda for all sites
#'
analyse_stan_res.all_binom_04 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda", "beta", "gamma", "sigma_beta", "sigma_gamma") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous all interactions results analyses:
  analyse_stan_res.misc_all_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.all_binom_04 <- c("post_param_plots", stan_res_mods.misc_all_binom)

#' Analyse all interactions binomial with alpha, betas and lambdas
#'
analyse_stan_res.all_binom_05 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda_bar", "beta", "lambda", "sigma_beta", "sigma_lambda") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous all interactions results analyses:
  analyse_stan_res.misc_all_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.all_binom_05 <- c("post_param_plots", stan_res_mods.misc_all_binom)

#' Analyse all interactions binomial with alpha, betas and single lambda
#'
analyse_stan_res.all_binom_06 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda", "beta", "sigma_beta") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous all interactions results analyses:
  analyse_stan_res.misc_all_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.all_binom_06 <- c("post_param_plots", stan_res_mods.misc_all_binom)

#' Analyse all interactions binomial with alpha, betas, gammas product and lambdas
#'
analyse_stan_res.all_binom_07 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda_bar", "beta", "zgamma", "lambda", "sigma_beta", 
      "sigma_gamma", "sigma_lambda") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous all interactions results analyses:
  analyse_stan_res.misc_all_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.all_binom_07 <- c("post_param_plots", stan_res_mods.misc_all_binom)

#' Analyse all interactions binomial with alpha, betas, gammas prod and single lambda
#'
analyse_stan_res.all_binom_08 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda", "beta", "zgamma", "sigma_beta", "sigma_gamma") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous all interactions results analyses:
  analyse_stan_res.misc_all_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.all_binom_08 <- c("post_param_plots", stan_res_mods.misc_all_binom)

#' Analyse all interactions binomial with two terms for bioclimatic suitability
#'
analyse_stan_res.all_binom_09 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda", "beta", "gamma", "sigma_beta", "sigma_gamma") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous all interactions results analyses:
  analyse_stan_res.misc_all_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules)
}
stan_res_mods.all_binom_09 <- c("post_param_plots", stan_res_mods.misc_all_binom)

#' Analyse all interactions binomial with single lambda for all sites (1 origin)
#'
analyse_stan_res.all_binom_14 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda", "mu", "beta", "gamma", "sigma_beta", "sigma_gamma") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous all interactions results analyses:
  analyse_stan_res.misc_all_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules, include_origin = 1)
}
stan_res_mods.all_binom_14 <- c("post_param_plots", stan_res_mods.misc_all_binom)

#' Analyse all interactions binomial with two terms for bioclimatic suitability (1 origin)
#'
analyse_stan_res.all_binom_19 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda", "mu", "beta", "gamma", "sigma_beta", "sigma_gamma") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous all interactions results analyses:
  analyse_stan_res.misc_all_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules, include_origin = 1)
}
stan_res_mods.all_binom_19 <- c("post_param_plots", stan_res_mods.misc_all_binom)

#' Analyse all interactions binomial with single lambda for all sites (2 origins)
#'
analyse_stan_res.all_binom_24 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda", "mu", "beta", "gamma", "sigma_beta", "sigma_gamma") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous all interactions results analyses:
  analyse_stan_res.misc_all_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules, include_origin = 2)
}
stan_res_mods.all_binom_24 <- c("post_param_plots", stan_res_mods.misc_all_binom)

#' Analyse all interactions binomial with two terms for bioclimatic suitability (2 origins)
#'
analyse_stan_res.all_binom_29 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder, prev_modules
) {
  # Path the previous analyses modules:
  prev_path <- file.path(res_folder, "prev_analyses.txt")
  
  # Plot posterior distribution and true value of parameters:
  if (!"post_param_plots" %in% prev_modules) {
    c("alpha", "lambda", "mu", "beta", "gamma", "sigma_beta", "sigma_gamma") %>%
      stan_analyses_plot_save_params_post(rstan_fit, ., res_folder)
    readr::write_lines("post_param_plots", prev_path, append = TRUE)
  }
  
  # Perform miscellaneous all interactions results analyses:
  analyse_stan_res.misc_all_binom(spec, data, start, cmdstan_fit, rstan_fit, 
                                  res_folder, prev_modules, include_origin = 2)
}
stan_res_mods.all_binom_29 <- c("post_param_plots", stan_res_mods.misc_all_binom)
