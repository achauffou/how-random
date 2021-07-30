# Generic script to compare and analyse results from all Stan models ===========
#' Generic function to compare Stan simulations
#'
compare_stan_res <-function(specs, res_folder, ...) {
  # Create results folder if it does not exist:
  comp_folder <- file.path(res_folder, "comparisons")
  dir.create(comp_folder, showWarnings = FALSE, recursive = TRUE)

  # Compare WAIC and LOOIC of models with identical data:
  compare_all_waics(specs, res_folder)
  compare_all_looics(specs, res_folder)

  # Return the current system time:
  as.character(Sys.time())
}

# Scripts to compare WAIC of models ============================================
#' Compare WAIC of all models (one comparison per models group with same data)
#'
compare_all_waics <- function(specs, res_folder) {
  # Get all analyses data generation models and arguments:
  comp_groups <- sapply(specs, function(x) x$data_comp_group)
  unique_comp_groups <- comp_groups[duplicated(comp_groups)] %>% unique()

  # Apply WAIC comparison to all groups of models with the same data:
  purrr::map(unique_comp_groups, function(x) {
    the_specs <- specs[which(comp_groups == x)]
    compare_similar_waics(the_specs, res_folder)
  })
}

#' Compare WAIC of models that have identical data
#'
compare_similar_waics <- function(specs, res_folder) {
  # Get hash name of the data generation specifications:
  hash_name <- specs[[1]]$data_comp_group

  # Read WAIC objects:
  waics <- lapply(specs, function(spec) {
    readRDS(file.path(res_folder, spec$name, "waic.rds"))
  })
  names(waics) <- purrr::map_chr(specs, ~ .$name)
  
  # Perform, plot and save comparison:
  comp <- loo::loo_compare(waics)
  suppressMessages({plot(comp) %>%
    ggsave(file.path(
      res_folder, "comparisons", paste0("waic_", hash_name, ".pdf")
    ), ., device = "pdf")})
  saveRDS(comp, file = file.path(
    res_folder, "comparisons", paste0("waic_", hash_name, ".rds")
  ))
}


# Scripts to compare LOOIC of models ============================================
#' Compare LOOIC of all models (one comparison per models group with same data)
#'
compare_all_looics <- function(specs, res_folder) {
  # Get all analyses data generation models and arguments:
  comp_groups <- sapply(specs, function(x) x$data_comp_group)
  unique_comp_groups <- comp_groups[duplicated(comp_groups)] %>% unique()
  
  # Apply LOOIC comparison to all groups of models with the same data:
  purrr::map(unique_comp_groups, function(x) {
    the_specs <- specs[which(comp_groups == x)]
    compare_similar_looics(the_specs, res_folder)
  })
}

#' Compare LOOIC of models that have identical data
#'
compare_similar_looics <- function(specs, res_folder) {
  # Get hash name of the data generation specifications:
  hash_name <- specs[[1]]$data_comp_group
  
  # Read LOOIC objects:
  looics <- lapply(specs, function(spec) {
    readRDS(file.path(res_folder, spec$name, "looic.rds"))
  })
  names(looics) <- purrr::map_chr(specs, ~ .$name)
  
  # Perform, plot and save comparison:
  comp <- loo::loo_compare(looics)
  suppressMessages({plot(comp) %>%
    ggsave(file.path(
      res_folder, "comparisons", paste0("looic_", hash_name, ".pdf")
    ), ., device = "pdf")})
  saveRDS(comp, file = file.path(
    res_folder, "comparisons", paste0("looic_", hash_name, ".rds")
  ))
}


# Plot models comparison =======================================================
#' Method to plot objects of class compare.loo
#' 
plot.compare.loo <- function(comp, keep = NULL) {
  # Get name of the information criterion:
  IC_name <- toupper(colnames(comp)[7])
  
  # If models to keep are specified, drop the others:
  if (!is.null(keep)) {
    comp <- comp[keep,]
  }
  
  # Prepare plot data:
  plt_data <- data.table::data.table(
    model = rownames(comp),
    value = comp[, 7],
    se_value = comp[, 8],
    diff = -2 * comp[, 1],
    se_diff = 2 * comp[, 2]
  )
  plt_data[1, ':='(diff = NA, se_diff = NA)]
  plt_data[, ':='(model_grp = -.GRP), by = .(model)]
  plt_data[, ':='(
    se_value_min = value - se_value,
    se_value_max = value + se_value,
    diff_value = min(value) + diff,
    se_diff_value_min = min(value) + diff - se_diff,
    se_diff_value_max = min(value) + diff + se_diff
  )]
  
  # Plot comparison:
  ggplot(plt_data) +
    geom_vline(aes(xintercept = min(plt_data$value)), color = "grey", size= 0.5) +
    geom_point(aes(x = value, y = model_grp), shape = 1) +
    geom_linerange(aes(xmin = se_value_min, xmax = se_value_max, y = model_grp)) +
    geom_point(
      aes(x = diff_value, y = model_grp + 0.5), 
      data = plt_data[-1], 
      shape = 2,
      color = "grey"
    ) +
    geom_linerange(
      aes(xmin = se_diff_value_min, xmax = se_diff_value_max, y = model_grp + 0.5), 
      data = plt_data[-1],
      color = "grey"
    ) +
    scale_x_continuous(name = "deviance") +
    scale_y_continuous(name = NULL, breaks = plt_data$model_grp, labels = plt_data$model) +
    ggtitle(IC_name) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5), 
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey", linetype = "dotted", size = 0.5),
      panel.grid.minor.y = element_blank()
    )
}
