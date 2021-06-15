# Overall functions to generate simulation data and starting values ============
#' Generate data for a STAN model simulation from its specification
#' 
generate_stan_sim_data <- function(spec, results_folder = "results/stan_sim") {
  # Create results folder if it does not exist:
  sim_folder <- file.path(results_folder, spec$name)
  dir.create(sim_folder, showWarnings = FALSE, recursive = TRUE)
  
  # If the previous run had the same specification, return its file:
  out_path <- file.path(sim_folder, "data.rds")
  last_spec_path <- file.path(sim_folder, "last_spec.rds")
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
  paste0("generate_stan_data.", spec$data_generation_model) %>%
    get() %>%
    do.call(args = fun_args) %>%
    saveRDS(file = out_path)
  out_path
}

#' Generate starting values for a STAN model simulation from its specification
#' 
generate_stan_sim_start_values <- function(spec, results_folder = "results/stan_sim") {
  # Create results folder if it does not exist:
  sim_folder <- file.path(results_folder, spec$name)
  dir.create(sim_folder, showWarnings = FALSE, recursive = TRUE)
  
  # If the previous run had the same specification, return its file:
  out_path <- file.path(sim_folder, "start_values.rds")
  last_spec_path <- file.path(sim_folder, "last_spec.rds")
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
  paste0("generate_stan_start_values.", spec$data_generation_model) %>%
    get() %>%
    do.call(args = fun_args) %>%
    saveRDS(file = out_path)
  out_path
}


# Custom distributions =========================================================
#' Draw random numbers from a truncated beta distribution
#' 
rbetacut <- function(x, shape1, shape2, ncp = 0, low_cut = 0, high_cut = 1) {
  single_rbeta_cut <- function(x, shape1, shape2, ncp, low_cut, high_cut) {
    out <- 1.2
    while (out < low_cut | out > high_cut) {
      out <- rbeta(1, shape1, shape2, ncp)
    }
    (out - low_cut) / (high_cut - low_cut)
  }
  sapply(1:x, single_rbeta_cut, shape1, shape2, ncp, low_cut, high_cut)
}


# Specific functions to generate data and starting values ======================
#' Generate data for pollination binomial with centered intercepts only
#' 
generate_stan_data.pol_binom_01 <- function(nb_sites, nb_pla, nb_pol) {
  # Create data.table with all interactions:
  data <- expand.grid(site_id = 1:nb_sites, pla_id = 1:nb_pla, pol_id = 1:nb_pol) %>%
    data.table::as.data.table()
  
  # Sample parameters (center intercept parameters):
  beta <- rnorm(nb_sites, 0, 1) %>% scale(scale = FALSE)
  gamma_pla <- rnorm(nb_pla, 0, 1) %>% scale(scale = FALSE)
  gamma_pol <- rnorm(nb_pol, 0, 1) %>% scale(scale = FALSE)
  
  # Compute response variable:
  data[, p := boot::inv.logit(
    beta[site_id] + gamma_pla[pla_id] + gamma_pol[pol_id]
  )]
  data[, Y := purrr::map_int(p, ~rbinom(1, 1, .))]
  
  # Return data specified as a list:
  list(
    nb_sites = nb_sites,
    nb_pla = nb_pla,
    nb_pol = nb_pol,
    nb_int = nrow(data),
    Y_array = data[, .(Y, site_id, pla_id, pol_id)],
    sigma_beta = sd(beta),
    sigma_gamma_pla = sd(gamma_pla),
    sigma_gamma_pol = sd(gamma_pol),
    beta = beta,
    gamma_pla = gamma_pla,
    gamma_pol = gamma_pol
  )
}

generate_stan_start_values.pol_binom_01 <- function(nb_sites, nb_pla, nb_pol) {
  list(
    zbeta = rep(0, nb_sites),
    zgamma_pla = rep(0, nb_pla),
    zgamma_pol = rep(0, nb_pol),
    sigma_beta = 0.1,
    sigma_gamma_pla = 0.1,
    sigma_gamma_pol = 0.1
  )
}

#' Generate data for pollination binomial with intercepts only
#' 
generate_stan_data.pol_binom_02 <- function(nb_sites, nb_pla, nb_pol) {
  # Create data.table with all interactions:
  data <- expand.grid(site_id = 1:nb_sites, pla_id = 1:nb_pla, pol_id = 1:nb_pol) %>%
    data.table::as.data.table()
  
  # Sample parameters (center intercept parameters):
  alpha <- rnorm(1, 0, 1)
  beta <- rnorm(nb_sites, 0, 1) %>% scale(scale = FALSE)
  gamma_pla <- rnorm(nb_pla, 0, 1) %>% scale(scale = FALSE)
  gamma_pol <- rnorm(nb_pol, 0, 1) %>% scale(scale = FALSE)
  
  # Compute response variable:
  data[, p := boot::inv.logit(
    alpha + beta[site_id] + gamma_pla[pla_id] + gamma_pol[pol_id]
  )]
  data[, Y := purrr::map_int(p, ~rbinom(1, 1, .))]
  
  # Return data specified as a list:
  list(
    nb_sites = nb_sites,
    nb_pla = nb_pla,
    nb_pol = nb_pol,
    nb_int = nrow(data),
    Y_array = data[, .(Y, site_id, pla_id, pol_id)],
    alpha = alpha,
    sigma_beta = sd(beta),
    sigma_gamma_pla = sd(gamma_pla),
    sigma_gamma_pol = sd(gamma_pol),
    beta = beta,
    gamma_pla = gamma_pla,
    gamma_pol = gamma_pol
  )
}

generate_stan_start_values.pol_binom_02 <- function(nb_sites, nb_pla, nb_pol) {
  list(
    alpha = 0,
    zbeta = rep(0, nb_sites),
    zgamma_pla = rep(0, nb_pla),
    zgamma_pol = rep(0, nb_pol),
    sigma_beta = 0.1,
    sigma_gamma_pla = 0.1,
    sigma_gamma_pol = 0.1
  )
}
