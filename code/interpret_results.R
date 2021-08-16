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
    "Non-neutral, two unpooled CS terms",
    "Non-neutral, one pooled CS term",
    "Non-neutral, no CS term",
    "Non-neutral, two pooled CS terms",
    "Non-neutral, one unpooled CS term",
    "Neutral, one pooled CS term",
    "Neutral, one unpooled CS term",
    "Neutral, two unpooled CS terms",
    "Neutral, no CS term",
    "Null model, per-type intercept only"
  )
  suppressMessages({suppressWarnings({
    plt1 <- plot(plt_data1)
    ggsave(file_path[1], plt1, device = "pdf", ...)
  })})
  
  # Non-neutral models without origin status:
  suppressMessages({suppressWarnings({
    plt2 <- plot(plt_data1, keep = 1:5)
    ggsave(file_path[2], plt2, device = "pdf", ...)
  })})
  
  # Non-neutral models with origin status:
  plt_data2 <- readRDS(file.path(comp_folder, "looic_all_bioclim_origin1.rds"))
  rownames(plt_data2) <- c(
    "No origin status, two unpooled CS terms",
    "Origin status, two unpooled CS terms",
    "No origin status, one unpooled CS term",
    "Origin status, one unpooled CS term",
    "Origin status, neutral, two unpooled CS terms"
  )
  suppressMessages({suppressWarnings({
    plt3 <- plot(plt_data2, keep = 1:4)
    ggsave(file_path[3], plt3, device = "pdf", ...)
  })})
  list(plt1, plt2, plt3)
}
