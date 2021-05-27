# Compile and run a Stan model =================================================
#' Compile, run and save the outcome in a given folder
#' 
run_stan_model <- function(
  spec, data, start_vals, src_folder, bin_folder, res_folder
) {
  # Create bin folder if it does not exist:
  dir.create(bin_folder, showWarnings = FALSE, recursive = TRUE)
  
  # Create results folder if it does not exist:
  out_folder <- file.path(res_folder, spec$name)
  dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
  
  # Compile the Stan model:
  model_path <- paste0(src_folder, "/", spec$stan_model, ".stan")
  model <- cmdstanr::cmdstan_model(
    model_path, cpp_options = list(stan_threads = TRUE), dir = bin_folder
  )
  
  # Prepare sampling arguments:
  chains <- spec$chains
  nb_cores <- parallel::detectCores()
  parallel_chains <- min(chains, nb_cores)
  threads_per_chain <- floor(nb_cores / chains)
  if (threads_per_chain < 1) threads_per_chain <- 1
  initial_values <- lapply(1:chains, function(x) start_vals)
  fun_args <- c(list(
    data = data,
    init = initial_values,
    chains = chains,
    threads_per_chain = threads_per_chain,
    parallel_chains = parallel_chains
  ), spec$data_sample_args)
  
  # Sample the compiled Stan model:
  fit <- do.call(model$sample, args = fun_args)
  
  # Save cmdStanMCMC and RStan fit objects to results folder:
  out_path1 <- file.path(out_folder, "cmdstan-fit.rds")
  saveRDS(fit, file = out_path1)
  out_path2 <- file.path(out_folder, "rstan-fit.rds")
  saveRDS(rstan::read_stan_csv(fit$output_files()), out_path2)
  c(cmdstan_fit = out_path1, rstan_fit = out_path2)
}