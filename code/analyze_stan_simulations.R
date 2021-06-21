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
stan_analyses_auc <- function(y, link, res_folder) {
  # Compute and save the Area Under Curve values for all posterior draws:
  if (!file.exists(file.path(res_folder, "auc.rds"))) {
    auc <- 1:dim(link)[1] %>%
      parallel::mclapply(function(x) {
        pROC::auc(y, link[x,], quiet = TRUE)
      }, mc.cores = get_nb_cpus()) %>%
      unlist()
    saveRDS(auc, file = file.path(res_folder, "auc.rds"))
  } else {
    auc <- readRDS(file.path(res_folder, "auc.rds"))
  }

  # Plot and save the distribution of the AUC:
  plt <- ggplot(data.frame(auc = auc), aes(x = auc)) +
    geom_density() +
    xlab("Area Under Curve") +
    ylab("Density")
    theme_classic()
  suppressMessages({ggsave(file.path(res_folder, "auc.pdf"), plt, device = "pdf")})
}

#' Plot the ROC curve of a Bayesian binomial model
#'
stan_analyses_roc <- function(y, link, res_folder) {
  # Compute and save the ROC curve coordinates for all posterior draws:
  if (!file.exists(file.path(res_folder, "roc.rds"))) {
    roc <- 1:dim(link)[1] %>%
      parallel::mclapply(function(x) {
        coords <- y %>%
          pROC::roc(link[x,], auc = FALSE, ci = FALSE, quiet = TRUE) %>%
          pROC::coords(x = seq(0, 1, by = 0.001), input = "specificity")
        coords$sample <- x
        coords
      }, mc.cores = get_nb_cpus()) %>% data.table::rbindlist()
    data.table::fwrite(roc, file.path(res_folder, "roc.csv"))
  } else {
    roc <- data.table::fread(file.path(res_folder, "roc.csv"))
  }

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
  suppressMessages({ggsave(file.path(res_folder, "roc.pdf"), plt, device = "pdf")})
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

  # Compute link:
  link <- link.pol_binom_01(data, rstan_fit)

  # Compute and plot AUC/ROC:
  stan_analyses_auc(data$Y_array$Y, link, res_folder)
  stan_analyses_roc(data$Y_array$Y, link, res_folder)
}

#' Link function of pollination binomial with centered intercepts only
#'
link.pol_binom_01 <- function(data, fit) {
  beta <- as.matrix(fit, pars = "beta")
  gamma_pla <- as.matrix(fit, pars = "gamma_pla")
  gamma_pol <- as.matrix(fit, pars = "gamma_pol")
  parallel::mclapply(1:(dim(fit)[1] * dim(fit)[2]), function(x) {
    p <- boot::inv.logit(
      beta[x, data$Y_array$site_id] + gamma_pla[x, data$Y_array$pla_id] +
      gamma_pol[x, data$Y_array$pol_id]
    )
    names(p) <- NULL
    p
  }, mc.cores = get_nb_cpus()) %>% do.call(rbind, args = .)
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

  # Compute link:
  link <- link.pol_binom_02(data, rstan_fit)

  # Compute and plot AUC/ROC:
  stan_analyses_auc(data$Y_array$Y, link, res_folder)
  stan_analyses_roc(data$Y_array$Y, link, res_folder)
}

#' Link function of pollination binomial with intercepts only
#'
link.pol_binom_02 <- function(data, fit) {
  alpha <- as.matrix(fit, pars = "alpha")
  beta <- as.matrix(fit, pars = "beta")
  gamma_pla <- as.matrix(fit, pars = "gamma_pla")
  gamma_pol <- as.matrix(fit, pars = "gamma_pol")
  parallel::mclapply(1:(dim(fit)[1] * dim(fit)[2]), function(x) {
    p <- boot::inv.logit(
      alpha[x] + beta[x, data$Y_array$site_id] +
        gamma_pla[x, data$Y_array$pla_id] + gamma_pol[x, data$Y_array$pol_id]
    )
    names(p) <- NULL
    p
  }, mc.cores = get_nb_cpus()) %>% do.call(rbind, args = .)
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

  # Compute link:
  link <- link.pol_binom_03(data, rstan_fit)

  # Compute and plot AUC/ROC:
  stan_analyses_auc(data$Y_array$Y, link, res_folder)
  stan_analyses_roc(data$Y_array$Y, link, res_folder)
}

#' Link function of pollination binomial with intercepts and slope
#'
link.pol_binom_03 <- function(data, fit) {
  alpha <- as.matrix(fit, pars = "alpha")
  beta <- as.matrix(fit, pars = "beta")
  gamma_pla <- as.matrix(fit, pars = "gamma_pla")
  gamma_pol <- as.matrix(fit, pars = "gamma_pol")
  lambda <- as.matrix(fit, pars = "lambda")
  parallel::mclapply(1:(dim(fit)[1] * dim(fit)[2]), function(x) {
    p <- boot::inv.logit(
      alpha[x] + beta[x, data$Y_array$site_id] +
        gamma_pla[x, data$Y_array$pla_id] + gamma_pol[x, data$Y_array$pol_id] +
        lambda[x, data$Y_array$site_id] * data$SS
    )
    names(p) <- NULL
    p
  }, mc.cores = get_nb_cpus()) %>% do.call(rbind, args = .)
}
