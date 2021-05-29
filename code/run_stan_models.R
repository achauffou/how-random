# Compile and run a Stan model =================================================
#' Compile, run and save the outcome in a given folder
#' 
run_stan_model <- function(
  spec, data, start_vals, src_file, bin_folder, res_folder
) {
  # Create bin folder if it does not exist:
  dir.create(bin_folder, showWarnings = FALSE, recursive = TRUE)
  
  # Create results folder if it does not exist:
  out_folder <- file.path(res_folder, spec$name)
  dir.create(out_folder, showWarnings = FALSE, recursive = TRUE)
  
  # Compile the Stan model:
  model <- cmdstanr::cmdstan_model(
    src_file, cpp_options = list(stan_threads = TRUE), dir = bin_folder
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
  ), spec$stan_sample_args)
  
  # Sample the compiled Stan model:
  fit <- do.call(model$sample, args = fun_args)
  
  # Print and save diagnostic to a file:
  print(fit$cmdstan_diagnose())
  sink(file = file.path(out_folder, "cmdstan-diagnostic.txt"))
  print(fit$cmdstan_diagnose())
  sink()
  
  # Save cmdStanMCMC and RStan fit objects to results folder:
  out_path1 <- file.path(out_folder, "cmdstanr-fit.rds")
  fit$save_object(file = out_path1)
  out_path2 <- file.path(out_folder, "rstan-fit.rds")
  saveRDS(rstan::read_stan_csv(fit$output_files()), out_path2)
  c(cmdstan_fit = out_path1, rstan_fit = out_path2)
}
