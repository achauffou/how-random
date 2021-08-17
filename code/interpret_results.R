# Functions that provide additional results interpretations ====================
#' Plot a graph illustrating the sensitivity analysis for climate suitability
#' 
plot_bioclim_sensitivity_illustration <- function(errors, file_path, case = "indiv", ...) {
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  plt_data <- data.table::rbindlist(lapply(1:5, function(x) errors[[x]][[case]]))
  suppressMessages({suppressWarnings({
    plt <- ggplot(plt_data, aes(x = log(sample_size), y = mae)) +
      geom_point() +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      geom_hline(yintercept = 0.15, color = "red", size = 0.5) +
      scale_x_continuous(name = "Log number of occurrences") +
      scale_y_continuous(name = "Average absolute error", limits = c(0, 0.6), expand = c(0,NA)) +
      theme_classic()
    ggsave(file_path, plt, device = "pdf", ...)
  })})
  plt
}

#' Plot graphs for the models comparison
#' 
plot_models_comp <- function(file_path, comp_folder, last_comp_update, ...) {
  dirname(file_path) %>% 
    purrr::map(~ dir.create(., showWarnings = FALSE, recursive = TRUE))
  
  # All models without origin status:
  plt_data1 <- readRDS(file.path(comp_folder, "looic_all_bioclim.rds"))
  rownames(plt_data1) <- c(
    "(5c) non-neutral, two unpooled CS terms",
    "(5b) non-neutral, one pooled CS term",
    "(3) non-neutral, no CS term",
    "(5d) non-neutral, two pooled CS terms",
    "(5a) non-neutral, one unpooled CS term",
    "(4b) neutral, one pooled CS term",
    "(4a) neutral, one unpooled CS term",
    "(4c) neutral, two unpooled CS terms",
    "(2) neutral, no CS term",
    "(1) null model, per-type intercept only"
  )
  suppressMessages({suppressWarnings({
    plt1 <- plot(plt_data1)
    ggsave(file_path[1], plt1, device = "pdf", ...)
  })})
  
  # Non-neutral models without origin status:
  suppressMessages({suppressWarnings({
    plt2 <- plot(plt_data1, keep = 1:5)
    ggsave(file_path[2], plt2, device = "pdf", width = 12, height = 6, units = "cm")
  })})
  
  # Non-neutral models with origin status:
  plt_data2 <- readRDS(file.path(comp_folder, "looic_all_bioclim_origin1.rds"))
  rownames(plt_data2) <- c(
    "(5c) no origin status, two unpooled CS terms",
    "(6c) origin status, two unpooled CS terms",
    "(5a) no origin status, one unpooled CS term",
    "(6a) origin status, one unpooled CS term",
    "(7c) origin status, neutral, two unpooled CS terms"
  )
  suppressMessages({suppressWarnings({
    plt3 <- plot(plt_data2, keep = 1:4)
    ggsave(file_path[3], plt3, device = "pdf", width = 12, height = 5, units = "cm")
  })})
  list(plt1, plt2, plt3)
}

#' Plot the posterior distribution of sigma_beta and sigma_gamma
#' 
plot_sigma_beta_gamma_distr <- function(fit_path, file_path, ...) {
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  fit <- readRDS(fit_path)
  suppressMessages({suppressWarnings({
    plt <- rstan::plot(fit, pars = c("sigma_beta", "sigma_gamma")) +
      scale_y_continuous(breaks = 1:6, expand = c(0.05, 0.05), labels = c(
        latex2exp::TeX("Seed dispersers $\\sigma_\\gamma$"),
        latex2exp::TeX("Pollinators $\\sigma_\\gamma$"),
        latex2exp::TeX("Plants (seed dispersal) $\\sigma_\\gamma$"),
        latex2exp::TeX("Plants (pollination) $\\sigma_\\gamma$"),
        latex2exp::TeX("Seed dispersal $\\sigma_\\beta$"),
        latex2exp::TeX("Pollination $\\sigma_\\beta$")
      )) +
      theme(plot.margin = grid::unit(c(5, 5, 5, 5), unit = "mm"))
    ggsave(file_path, plt, device = "pdf", width = 12, height = 6, units = "cm")
  })})
  file_path
}

#' Plot the posterior distribution of lambda parameters
#' 
plot_lambda_distr <- function(fit_path, file_path, ...) {
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  fit <- readRDS(fit_path)
  suppressMessages({suppressWarnings({
    plt <- rstan::plot(fit, pars = c("lambda")) +
      scale_y_continuous(breaks = 1:4, expand = c(0.05, 0.05), labels = c(
        latex2exp::TeX("Seed dispersers $\\lambda$"),
        latex2exp::TeX("Pollinators $\\lambda$"),
        latex2exp::TeX("Plants (seed dispersal) $\\lambda$"),
        latex2exp::TeX("Plants (pollination) $\\lambda$")
      )) +
      theme(plot.margin = grid::unit(c(5, 5, 5, 5), unit = "mm"))
    ggsave(file_path, plt, device = "pdf", width = 12, height = 6, units = "cm")
  })})
  file_path
}

#' Plot the posterior distribution of eta parameters
#' 
plot_eta_distr <- function(fit_path, file_path, ...) {
  dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)
  fit <- readRDS(fit_path)
  suppressMessages({suppressWarnings({
    plt <- rstan::plot(fit, pars = c("mu")) +
      scale_y_continuous(breaks = 1:2, expand = c(0.1, 0.1), labels = c(
        latex2exp::TeX("Plants (seed dispersal) $\\eta$"),
        latex2exp::TeX("Plants (pollination) $\\eta$")
      )) +
      theme(plot.margin = grid::unit(c(5, 5, 5, 5), unit = "mm"))
    ggsave(file_path, plt, device = "pdf", width = 12, height = 3, units = "cm")
  })})
  file_path
}
