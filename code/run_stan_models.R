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
  
  # Last specification file and output files:
  out_path1 <- file.path(out_folder, "cmdstanr-fit.rds")
  out_path2 <- file.path(out_folder, "rstan-fit.rds")
  last_spec_path <- file.path(out_folder, "last_spec.rds")
  last_src_md5sum_path <- file.path(out_folder, "last_src_md5sum.txt")
  
  # If the previous run had the same specification, return its files:
  if (
    file.exists(last_spec_path) & file.exists(out_path1) & 
    file.exists(out_path2) & file.exists(last_src_md5sum_path)
  ) {
    if (identical(readRDS(last_spec_path), spec)) {
      con <- file(last_src_md5sum_path, "r")
      last_md5sum <- readLines(con, n=1)
      close(con)
      this_md5sum <- tools::md5sum(src_file)
      if (last_md5sum == this_md5sum) {
        message(paste(spec$name, "results seem already up-to-date,", 
                      "skipping HMC sampling."))
        return(c(cmdstan_fit = out_path1, rstan_fit = out_path2))
      }
    }
  }
  
  # Compile the Stan model:
  model <- cmdstanr::cmdstan_model(
    src_file, cpp_options = list(stan_threads = TRUE), dir = bin_folder
  )
  
  # Prepare sampling arguments:
  chains <- spec$chains
  nb_cores <- get_nb_cpus()
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
  
  # Remove last_analysis.txt if it exists to ensure that analysis is performed:
  if (file.exists(file.path(out_folder, "last_analysis.txt"))) {
    file.remove(file.path(out_folder, "last_analysis.txt"))
  }
  
  # Print and save diagnostic to a file:
  print(fit$cmdstan_diagnose())
  sink(file = file.path(out_folder, "cmdstan-diagnostic.txt"))
  print(fit$cmdstan_diagnose())
  sink()
  
  # Save cmdStanMCMC and RStan fit objects to results folder:
  fit$save_object(file = out_path1)
  saveRDS(rstan::read_stan_csv(fit$output_files()), out_path2)
  
  # Save specification and source md5sum to skip this step in the future:
  saveRDS(spec, last_spec_path)
  con <- file(last_src_md5sum_path, "w")
  writeLines(tools::md5sum(src_file), con)
  close(con)
  
  # Return file names of CmdStanMCMC and Stanfit objects:
  c(cmdstan_fit = out_path1, rstan_fit = out_path2)
}
