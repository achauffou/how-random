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
#' Generate data for a pollination logit model with no alpha intercept
#' 
generate_stan_data.pol_logit_f <- function(nb_sites, nb_pla, nb_pol) {
  # Generate degree and optimal suitability:
  D_pla <- rbetacut(nb_pla, 3, 2, high_cut = 0.9)
  D_pol <- rbetacut(nb_pol, 3, 2, high_cut = 0.9)
  S_opt_pla <- runif(nb_pla, 0, 1)
  S_opt_pol <- runif(nb_pol, 0, 1)
  env_sit <- runif(nb_sites, 0, 1)
  S_pla <- 1 - abs(outer(S_opt_pla, env_sit, "-"))
  S_pol <- 1 - abs(outer(S_opt_pol, env_sit, "-"))
  
  # Create data.table with all interactions:
  data <- lapply(1:nb_sites, function(x) {
    prop_sample <- rbeta(1, 2, 6)
    data <- expand.grid(
      pla_id = sample(1:nb_pla, round(nb_pla * prop_sample)), 
      pol_id = sample(1:nb_pol, round(nb_pol * prop_sample)), 
      site_id = x
    ) %>% data.table::as.data.table()
  }) %>% data.table::rbindlist()
  data[, ':='(
    D = D_pla[pla_id] * D_pol[pol_id],
    S = S_pla[cbind(pla_id, site_id)] * S_pol[cbind(pol_id, site_id)]
  )]
  
  # Sample parameters:
  beta <- rnorm(nb_sites, 0, 1)
  gamma_pla <- rnorm(nb_pla, 0, 1)
  gamma_pol <- rnorm(nb_pol, 0, 1)
  lambda <- rnorm(nb_sites, 0, 1)
  nu <- rnorm(nb_sites, 0, 1)
  
  # Compute response variable:
  data[, p := boot::inv.logit(
    beta[site_id] + gamma_pla[pla_id] + gamma_pol[pol_id] + 
      lambda[site_id] * S + nu[site_id] * D
  )]
  data[, Y := purrr::map_int(p, ~rbinom(1, 1, .))]
  
  # Return data specified as a list:
  list(
    nb_sites = nb_sites,
    nb_pla = nb_pla,
    nb_pol = nb_pol,
    nb_int = nrow(data),
    Y_array = data[, .(Y, site_id, pla_id, pol_id)],
    S_pla = S_pla,
    S_pol = S_pol,
    D_pla = D_pla,
    D_pol = D_pol,
    beta_site = beta,
    gamma_pla = gamma_pla,
    gamma_pol = gamma_pol,
    lambda = lambda,
    nu = nu
  )
}

generate_stan_start_values.pol_logit_f <- function(nb_sites, nb_pla, nb_pol) {
  list(
    zbeta = rep(0, nb_sites),
    zgamma_pol = rep(0, nb_pol),
    zgamma_pla = rep(0, nb_pla),
    zlambda = rep(0, nb_sites),
    znu = rep(0, nb_sites),
    sigma_beta = 0.1,
    sigma_gamma_pla = 0.1,
    sigma_gamma_pol = 0.1,
    sigma_lambda = 0.1,
    sigma_nu = 0.1
  )
}

#' Generate data for a pollination logit model with a random alpha intercept
#' 
generate_stan_data.pol_logit_g <- function(nb_sites, nb_pla, nb_pol) {
  # Generate degree and optimal suitability:
  D_pla <- rbetacut(nb_pla, 3, 2, high_cut = 0.9)
  D_pol <- rbetacut(nb_pol, 3, 2, high_cut = 0.9)
  S_opt_pla <- runif(nb_pla, 0, 1)
  S_opt_pol <- runif(nb_pol, 0, 1)
  env_sit <- runif(nb_sites, 0, 1)
  S_pla <- 1 - abs(outer(S_opt_pla, env_sit, "-"))
  S_pol <- 1 - abs(outer(S_opt_pol, env_sit, "-"))
  
  # Create data.table with all interactions:
  data <- lapply(1:nb_sites, function(x) {
    prop_sample <- rbeta(1, 2, 6)
    data <- expand.grid(
      pla_id = sample(1:nb_pla, round(nb_pla * prop_sample)), 
      pol_id = sample(1:nb_pol, round(nb_pol * prop_sample)), 
      site_id = x
    ) %>% data.table::as.data.table()
  }) %>% data.table::rbindlist()
  data[, ':='(
    D = D_pla[pla_id] * D_pol[pol_id],
    S = S_pla[cbind(pla_id, site_id)] * S_pol[cbind(pol_id, site_id)]
  )]
  
  # Sample parameters:
  alpha <- rnorm(1, 0, 1)
  beta <- rnorm(nb_sites, 0, 1)
  gamma_pla <- rnorm(nb_pla, 0, 1)
  gamma_pol <- rnorm(nb_pol, 0, 1)
  lambda <- rnorm(nb_sites, 0, 1)
  nu <- rnorm(nb_sites, 0, 1)
  
  # Compute response variable:
  data[, p := boot::inv.logit(
    alpha + beta[site_id] + gamma_pla[pla_id] + gamma_pol[pol_id] + 
      lambda[site_id] * S + nu[site_id] * D
  )]
  data[, Y := purrr::map_int(p, ~rbinom(1, 1, .))]
  
  # Return data specified as a list:
  list(
    nb_sites = nb_sites,
    nb_pla = nb_pla,
    nb_pol = nb_pol,
    nb_int = nrow(data),
    Y_array = data[, .(Y, site_id, pla_id, pol_id)],
    S_pla = S_pla,
    S_pol = S_pol,
    D_pla = D_pla,
    D_pol = D_pol,
    alpha = alpha,
    beta_site = beta,
    gamma_pla = gamma_pla,
    gamma_pol = gamma_pol,
    lambda = lambda,
    nu = nu
  )
}

generate_stan_start_values.pol_logit_g <- function(nb_sites, nb_pla, nb_pol) {
  list(
    alpha = 0,
    zbeta = rep(0, nb_sites),
    zgamma_pol = rep(0, nb_pol),
    zgamma_pla = rep(0, nb_pla),
    zlambda = rep(0, nb_sites),
    znu = rep(0, nb_sites),
    sigma_beta = 0.1,
    sigma_gamma_pla = 0.1,
    sigma_gamma_pol = 0.1,
    sigma_lambda = 0.1,
    sigma_nu = 0.1
  )
}
