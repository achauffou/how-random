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
  if (!file.exists(file.path(res_folder, "last_analysis.txt"))) {
    fun_name <- paste0("analyse_stan_sim.", spec$stan_model)
    if (exists(fun_name)) {
      fun_name %>%
        get() %>%
        do.call(args = list(spec, data, start, cmdstan_fit, rstan_fit, res_folder))
      con <- file(file.path(res_folder, "last_analysis.txt"))
      writeLines(as.character(Sys.time()), con)
      close(con)
    }
  } else {
    message(paste("Results of simulation", spec$name, 
                  "seem to already have been analysed.",
                  "Skipping Stan results analysis..."))
  }

  # Return the current system time:
  as.character(Sys.time())
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
      ggsave(file.path(res_folder, paste0("param_post_true_", x, ".pdf")), .,
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
      ggsave(file.path(res_folder, paste0("param_error_", params[x], ".pdf")), .,
             device = "pdf")
  })})})
}

#' Compute and plot the AUC of a Bayesian binomial model
#'
stan_analyses_auc <- function(
  y, link, res_folder, auc_name = "auc", nb_samples = NULL
) {
  # Sample link:
  if(!is.null(nb_samples)) {
    link <- link[sample(1:dim(link)[1], nb_samples), ]
  }
  
  # Compute and save the Area Under Curve values for all posterior draws:
  auc <- 1:dim(link)[1] %>%
    parallel::mclapply(function(x) {
      pROC::auc(y, link[x,], quiet = TRUE)
    }, mc.cores = get_nb_cpus()) %>%
    unlist()
  saveRDS(auc, file = file.path(res_folder, paste0(auc_name, ".rds")))

  # Plot and save the distribution of the AUC:
  plt <- ggplot(data.frame(auc = auc), aes(x = auc)) +
    geom_density() +
    xlab("Area Under Curve") +
    ylab("Density") +
    theme_classic()
  suppressMessages({ggsave(
    file.path(res_folder, paste0(auc_name, ".pdf")), plt, device = "pdf"
  )})
}

#' Plot the ROC curve of a Bayesian binomial model
#'
stan_analyses_roc <- function(
  y, link, res_folder, roc_name = "roc", nb_samples = NULL
) {
  # Sample link:
  if(!is.null(nb_samples)) {
    link <- link[sample(1:dim(link)[1], nb_samples), ]
  }
  
  # Compute and save the ROC curve coordinates for all posterior draws:
  roc <- 1:dim(link)[1] %>%
    parallel::mclapply(function(x) {
      coords <- y %>%
        pROC::roc(link[x,], auc = FALSE, ci = FALSE, quiet = TRUE) %>%
        pROC::coords(x = seq(0, 1, by = 0.001), input = "specificity")
      coords$sample <- x
      coords
    }, mc.cores = get_nb_cpus()) %>% data.table::rbindlist()
  data.table::fwrite(roc, file.path(res_folder, paste0(roc_name, ".csv")))

  # Plot and save the distribution of the AUC:
  plt_data <- roc[, .(
    mean = mean(sensitivity),
    q05 = quantile(sensitivity, 0.05),
    q95 = quantile(sensitivity, 0.95)
  ), by = .(specificity)]
  plt <- ggplot(plt_data, aes(x = specificity, y = mean, ymin = q05, ymax = q95)) +
    geom_line() +
    geom_ribbon() +
    scale_x_reverse(name = "Specificity", limits = c(1, 0), expand = c(0, NA)) +
    scale_y_continuous(name = "Sensitivity", limits = c(0, 1), expand = c(0, NA)) +
    theme_classic()
  suppressMessages({ggsave(
    file.path(res_folder, paste0(roc_name, ".pdf")), plt, device = "pdf"
  )})
}


# Functions to analyse specific Stan simulations for pollination ===============
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
  c("alpha", "lambda_bar", "beta", "gamma_pla", "gamma_pol",
    "lambda", "sigma_beta", "sigma_gamma_pla", "sigma_gamma_pol",
    "sigma_lambda") %>%
    stan_sim_analyses_plot_save_params_post_true(rstan_fit, data, ., res_folder)

  # Plot posterior error of multilevel parameters:
  stan_sim_analyses_plot_save_params_errors(
    rstan_fit, data, c("beta", "gamma_pla", "gamma_pol", "lambda"),
    c("site_id", "pla_id", "pol_id", "site_id", "site_id"), res_folder
  )
}

#' Analyse pollination binomial with intercepts and seperate suitability terms
#'
analyse_stan_sim.pol_binom_09 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda_pla", "lambda_pol", "beta", "gamma_pla", "gamma_pol",
    "sigma_beta", "sigma_gamma_pla", "sigma_gamma_pol") %>%
    stan_sim_analyses_plot_save_params_post_true(rstan_fit, data, ., res_folder)
  
  # Plot posterior error of multilevel parameters:
  stan_sim_analyses_plot_save_params_errors(
    rstan_fit, data, c("beta", "gamma_pla", "gamma_pol"),
    c("site_id", "pla_id", "pol_id", "site_id", "site_id"), res_folder
  )
}

#' Analyse pollination binomial with intercepts and slopes
#'
analyse_stan_sim.pol_binom_24 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda", "mu_pla", "mu_pol", "beta", "gamma_pla", "gamma_pol",
    "sigma_beta", "sigma_gamma_pla", "sigma_gamma_pol") %>%
    stan_sim_analyses_plot_save_params_post_true(rstan_fit, data, ., res_folder)
  
  # Plot posterior error of multilevel parameters:
  stan_sim_analyses_plot_save_params_errors(
    rstan_fit, data, c("beta", "gamma_pla", "gamma_pol"),
    c("site_id", "pla_id", "pol_id"), res_folder
  )
}


# Functions to analyse specific Stan simulations for pollination ===============
#' Analyse all interactions binomial with intercepts and slopes
#'
analyse_stan_sim.all_binom_03 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda_bar", "beta", "gamma", "lambda", "sigma_beta",
    "sigma_gamma", "sigma_lambda") %>%
    stan_sim_analyses_plot_save_params_post_true(rstan_fit, data, ., res_folder)

  # Plot posterior error of multilevel parameters:
  stan_sim_analyses_plot_save_params_errors(
    rstan_fit, data, c("beta", "lambda"), c("site_id", "site_id"), res_folder
  )
  suppressWarnings({suppressMessages({
    nb_data <- rbind(
      data$Y_array[, .N, by = .(sp1_id)][, .(sp_id = sp1_id, N)],
      data$Y_array[, .N, by = .(sp2_id)][, .(sp_id = sp2_id, N)]
    )[order(sp_id)][['N']] %>%
      stan_sim_analyses_plot_param_error(rstan_fit, "gamma", data$gamma, .) %>%
      ggsave(file.path(res_folder, "param_error_gamma.pdf"), ., device = "pdf")
  })})
}

#' Analyse all interactions binomial with intercepts and single lambda slope
#'
analyse_stan_sim.all_binom_04 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda", "beta", "gamma", "sigma_beta", "sigma_gamma") %>%
    stan_sim_analyses_plot_save_params_post_true(rstan_fit, data, ., res_folder)
  
  # Plot posterior error of multilevel parameters:
  stan_sim_analyses_plot_save_params_errors(
    rstan_fit, data, c("beta"), c("site_id"), res_folder
  )
  suppressWarnings({suppressMessages({
    nb_data <- rbind(
      data$Y_array[, .N, by = .(sp1_id)][, .(sp_id = sp1_id, N)],
      data$Y_array[, .N, by = .(sp2_id)][, .(sp_id = sp2_id, N)]
    )[order(sp_id)][['N']] %>%
      stan_sim_analyses_plot_param_error(rstan_fit, "gamma", data$gamma, .) %>%
      ggsave(file.path(res_folder, "param_error_gamma.pdf"), ., device = "pdf")
  })})
}

#' Analyse all interactions binomial with intercepts and two lambda terms
#'
analyse_stan_sim.all_binom_09 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda", "beta", "gamma", "sigma_beta", "sigma_gamma") %>%
    stan_sim_analyses_plot_save_params_post_true(rstan_fit, data, ., res_folder)
  
  # Plot posterior error of multilevel parameters:
  stan_sim_analyses_plot_save_params_errors(
    rstan_fit, data, c("beta"), c("site_id"), res_folder
  )
  suppressWarnings({suppressMessages({
    nb_data <- rbind(
      data$Y_array[, .N, by = .(sp1_id)][, .(sp_id = sp1_id, N)],
      data$Y_array[, .N, by = .(sp2_id)][, .(sp_id = sp2_id, N)]
    )[order(sp_id)][['N']] %>%
      stan_sim_analyses_plot_param_error(rstan_fit, "gamma", data$gamma, .) %>%
      ggsave(file.path(res_folder, "param_error_gamma.pdf"), ., device = "pdf")
  })})
}

#' Analyse all interactions binomial with origin status
#'
analyse_stan_sim.all_binom_24 <- function(
  spec, data, start, cmdstan_fit, rstan_fit, res_folder
) {
  # Plot posterior distribution and true value of parameters:
  c("alpha", "lambda", "mu", "beta", "gamma", "sigma_beta", "sigma_gamma") %>%
    stan_sim_analyses_plot_save_params_post_true(rstan_fit, data, ., res_folder)
  
  # Plot posterior error of multilevel parameters:
  stan_sim_analyses_plot_save_params_errors(
    rstan_fit, data, c("beta"), c("site_id"), res_folder
  )
  suppressWarnings({suppressMessages({
    nb_data <- rbind(
      data$Y_array[, .N, by = .(sp1_id)][, .(sp_id = sp1_id, N)],
      data$Y_array[, .N, by = .(sp2_id)][, .(sp_id = sp2_id, N)]
    )[order(sp_id)][['N']] %>%
      stan_sim_analyses_plot_param_error(rstan_fit, "gamma", data$gamma, .) %>%
      ggsave(file.path(res_folder, "param_error_gamma.pdf"), ., device = "pdf")
  })})
}
