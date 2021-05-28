# Overall functions to generate simulation data and starting values ============
#' Generate data for a STAN model simulation from its specification
#' 
generate_stan_sim_data <- function(spec, results_folder = "results/stan_sim") {
  # Create results folder if it does not exist:
  sim_folder <- file.path(results_folder, spec$name)
  dir.create(sim_folder, showWarnings = FALSE, recursive = TRUE)
  
  # Generate and save data to the results folder:
  out_path <- file.path(sim_folder, "data.rds")
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
  
  # Generate and save starting values to the results folder:
  out_path <- file.path(sim_folder, "start_values.rds")
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
#' Simple model to generate plant-pollinator interaction data
#' 
generate_stan_data.pol_logit_a <- function(nb_pla, nb_pol, nb_sites) {
  # Generate degree and optimal suitability:
  K_pla <- rbetacut(nb_pla, 3, 2, high_cut = 0.9)
  K_pol <- rbetacut(nb_pol, 3, 2, high_cut = 0.9)
  S_opt_pla <- rnorm(nb_pla, 0, 1)
  S_opt_pol <- rnorm(nb_pol, 0, 1)
  env_sit <- rnorm(nb_sites, 0, 1)
  
  # Create data.table with all interactions:
  data <- expand.grid(
    pla_id = 1:nb_pla, pol_id = 1:nb_pol, site_id = 1:nb_sites) %>%
    data.table::as.data.table()
  data[, ':='(
    D = K_pla[pla_id] * K_pol[pol_id],
    S = abs(S_opt_pla[pla_id] - env_sit[site_id]) * 
      abs(S_opt_pol[pol_id] - env_sit[site_id])
  )]
  
  # Sample parameters:
  alpha <- rnorm(nb_sites, 0, 1)
  beta_pla <- rnorm(nb_pla, 0, 1)
  beta_pol <- rnorm(nb_pol, 0, 1)
  lambda <- rnorm(nb_sites, 0, 1)
  nu <- rnorm(nb_sites, 0, 1)
  
  # Compute response variable:
  data[, p := boot::inv.logit(
    alpha[site_id] + beta_pla[pla_id] + beta_pol[pol_id] + lambda[site_id] * D + 
      nu[site_id] * S
    )]
  data[, Y := purrr::map_int(p, ~rbinom(1, 1, .))]
  
  # Return data specified as a list:
  list(
    nb_pla = nb_pla,
    nb_pol = nb_pol,
    nb_sites = nb_sites,
    nb_int = nrow(data),
    Y = data[['Y']],
    site_id = data[['site_id']],
    pla_id = data[['pla_id']],
    pol_id = data[['pol_id']],
    D = data[['D']],
    S = data[['S']],
    alpha = alpha,
    beta_pol = beta_pol,
    beta_pla = beta_pla,
    lambda = lambda,
    nu = nu
  )
}

generate_stan_start_values.pol_logit_a <- function(nb_pla, nb_pol, nb_sites) {
  list(
    zalpha = rep(0, nb_sites),
    zbeta_pol = rep(0, nb_pol),
    zbeta_pla = rep(0, nb_pla),
    znu = rep(0, nb_sites),
    zlambda = rep(0, nb_sites),
    alpha_bar = 0,
    beta_pol_bar = 0,
    beta_pla_bar = 0,
    nu_bar = 0,
    lambda_bar = 0,
    sigma_a = 0.1,
    sigma_bp = 0.1,
    sigma_bpl = 0.1,
    sigma_n = 0.1,
    sigma_l = 0.1
  )
}
