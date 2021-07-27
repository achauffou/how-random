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
  
  # Perform and save comparison:
  loo::loo_compare(waics) %>%
    saveRDS(file = file.path(
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
  
  # Perform and save comparison:
  loo::loo_compare(looics) %>%
    saveRDS(file = file.path(
      res_folder, "comparisons", paste0("looic_", hash_name, ".rds")
    ))
}
