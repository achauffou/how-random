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
