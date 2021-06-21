# Overall functions to prepare Stan data and starting values ===================
#' Prepare data for a Stan model from its specification
#' 
prepare_stan_data <- function(spec, results_folder = "results/analyses") {
  # Create results folder if it does not exist:
  out_folder <- file.path(results_folder, spec$name)
  dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
  
  # If the previous run had the same specification, return its file:
  out_path <- file.path(out_folder, "data.rds")
  last_spec_path <- file.path(out_folder, "last_spec.rds")
  if (file.exists(last_spec_path) & file.exists(out_path)) {
    if (identical(readRDS(last_spec_path), spec)) {
      message(paste(spec$name, "results seem already up-to-date,", 
                    "skipping data generation."))
      return(out_path)
    }
  }
  
  # Generate and save data to the results folder:
  if (is.null(spec$data_generation_args)) {
    fun_args <- list()
  } else {
    fun_args <- spec$data_generation_args
  }
  paste0("prepare_stan_data.", spec$data_generation_model) %>%
    get() %>%
    do.call(args = c(interactions, fun_args)) %>%
    saveRDS(file = out_path)
  out_path
}

#' Generate starting values for a Stan model from its specification
#' 
prepare_stan_start_values <- function(spec, results_folder = "results/analyses") {
  # Create results folder if it does not exist:
  out_folder <- file.path(results_folder, spec$name)
  dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
  
  # If the previous run had the same specification, return its file:
  out_path <- file.path(out_folder, "start_values.rds")
  last_spec_path <- file.path(out_folder, "last_spec.rds")
  if (file.exists(last_spec_path) & file.exists(out_path)) {
    if (identical(readRDS(last_spec_path), spec)) {
      message(paste(spec$name, "results seem already up-to-date,", 
                    "skipping starting values generation."))
      return(out_path)
    }
  }
  
  # Generate and save starting values to the results folder:
  if (is.null(spec$data_generation_args)) {
    fun_args <- list()
  } else {
    fun_args <- spec$data_generation_args
  }
  paste0("prepare_stan_start_values.", spec$data_generation_model) %>%
    get() %>%
    do.call(args = fun_args) %>%
    saveRDS(file = out_path)
  out_path
}
