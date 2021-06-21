# Generic function to analyse Stan results =====================================
#' Determine and use the correct function to analyse a Stan outcome
#'
analyse_stan_res <- function(spec, data, start, fits) {
  # Determine results folder and read RDS files:
  res_folder <- dirname(fits[[1]])
  cmdstan_fit <- readRDS(fits[[1]])
  rstan_fit <- readRDS(fits[[2]])
  data <- readRDS(data)
  start <-  readRDS(start)
  
  # Analyse simulation outcome:
  fun_name <- paste0("analyse_stan_res.", spec$stan_model)
  if (exists(fun_name)) {
    fun_name %>%
      get() %>%
      do.call(args = list(spec, data, start, cmdstan_fit, rstan_fit, res_folder))
  }
  
  # Return RStan fit summary:
  list(rstan::summary(rstan_fit)$summary)
}
