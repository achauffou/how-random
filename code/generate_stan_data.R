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

#' Generate data for pollination binomial with intercepts and slopes
#' 
generate_stan_data.pol_binom_03 <- function(
  nb_sites, nb_pla, nb_pol, p_sample =  "1.0", rm_empty = TRUE
) {
  # Generate optimal suitability:
  S_opt_pla <- runif(nb_pla, 0, 1)
  S_opt_pol <- runif(nb_pol, 0, 1)
  env_sit <- runif(nb_sites, 0, 1)
  S_pla <- 1 - abs(outer(S_opt_pla, env_sit, "-"))
  S_pol <- 1 - abs(outer(S_opt_pol, env_sit, "-"))
  
  # Create data.table with all interactions:
  data <- lapply(1:nb_sites, function(x) {
    site_prop <- eval(parse(text = as.character(p_sample)))
    data <- expand.grid( 
      site_id = x,
      pla_id = sample(1:nb_pla, max(round(nb_pla * site_prop), 1)), 
      pol_id = sample(1:nb_pol, max(round(nb_pol * site_prop), 1))
    ) %>% data.table::as.data.table()
  }) %>% data.table::rbindlist()
  data[, ':='(
    S_pla = S_pla[cbind(pla_id, site_id)],
    S_pol = S_pol[cbind(pol_id, site_id)]
  )]
  
  # Standardize predictor variables:
  data[, ':='(
    SS = scale(S_pla * S_pol)
  )]
  
  # Sample parameters (center intercept parameters):
  alpha <- rbeta(1, 2, 4) * -1
  beta <- rnorm(nb_sites, 0, 1) %>% scale(scale = FALSE)
  gamma_pla <- rnorm(nb_pla, 0, 1) %>% scale(scale = FALSE)
  gamma_pol <- rnorm(nb_pol, 0, 1) %>% scale(scale = FALSE)
  lambda_bar <- rbeta(1, 2, 4)
  lambda <- rnorm(nb_sites, 0, 1) %>% scale(scale = FALSE) %>% add(lambda_bar)
  
  # Compute response variable:
  data[, p := boot::inv.logit(
    alpha + beta[site_id] + gamma_pla[pla_id] + gamma_pol[pol_id] +
      lambda[site_id] * SS
  )]
  data[, Y := purrr::map_int(p, ~rbinom(1, 1, .))]
  
  # Remove empty rows/columns:
  if (rm_empty == TRUE) {
    data[, sum_pla_ints := sum(Y), by = .(site_id, pla_id)]
    data[, sum_pol_ints := sum(Y), by = .(site_id, pol_id)]
    data <- data[sum_pla_ints > 0 & sum_pol_ints > 0]
  }
  
  # Return data specified as a list:
  list(
    nb_sites = nb_sites,
    nb_pla = nb_pla,
    nb_pol = nb_pol,
    nb_int = nrow(data),
    Y_array = data[, .(Y, site_id, pla_id, pol_id)],
    SS = data$SS,
    alpha = alpha,
    lambda_bar = lambda_bar,
    sigma_beta = sd(beta),
    sigma_gamma_pla = sd(gamma_pla),
    sigma_gamma_pol = sd(gamma_pol),
    sigma_lambda = sd(lambda),
    beta = beta,
    gamma_pla = gamma_pla,
    gamma_pol = gamma_pol,
    lambda = lambda,
    SS_mean = mean(data$S_pla * data$S_pol),
    SS_sd = sd(data$S_pla * data$S_pol)
  )
}

generate_stan_start_values.pol_binom_03 <- function(
  nb_sites, nb_pla, nb_pol, p_sample = "1.0", rm_empty = TRUE
) {
  list(
    alpha = 0,
    lambda_bar = 0,
    zbeta = rep(0, nb_sites),
    zgamma_pla = rep(0, nb_pla),
    zgamma_pol = rep(0, nb_pol),
    zlambda = rep(0, nb_sites),
    sigma_beta = 0.1,
    sigma_gamma_pla = 0.1,
    sigma_gamma_pol = 0.1,
    sigma_lambda = 0.1
  )
}

#' Generate data for all interaction types binomial with intercepts and slopes
#' 
generate_stan_data.all_binom_03 <- function(
  nb_sites, nb_spp, nb_types, p_sample =  "1.0", rm_empty = TRUE
) {
  # Assign all species and sites to an interaction type and group:
  site_type <- sample(1:nb_types, nb_sites, replace = TRUE)
  sp_group <- sample(1:(2 * nb_types), nb_spp, replace = TRUE)
  
  # Sample parameters:
  alpha <- rbeta(1, 4, 2) * -2
  lambda_bar <- rbeta(nb_types, 4, 2)
  sigma_beta <- rbeta(1, 4, 2) * 1.2
  sigma_gamma <- rbeta(2 * nb_types, 4, 2) * 1.2
  sigma_lambda <- rbeta(nb_types, 4, 2) * 1.2
  zbeta <- rnorm(nb_sites, 0, 1) %>% scale()
  zgamma <- rnorm(nb_spp, 0, 1) %>% scale()
  zlambda <- rnorm(nb_sites, 0, 1) %>% scale()
  
  # Compute non-centered parametrization:
  beta <- zbeta * sigma_beta
  gamma <- zgamma * sigma_gamma[sp_group]
  lambda <- zlambda * sigma_lambda[site_type] + lambda_bar[site_type]
  
  # Generate optimal suitability:
  S_opt <- runif(nb_spp, 0, 1)
  env_sit <- runif(nb_sites, 0, 1)
  S <- 1 - abs(outer(S_opt, env_sit, "-"))
  
  # Generate data interaction type by interaction type:
  data <- lapply(1:nb_types, function(type) {
    # Get IDs of species and sites that belong to this interaction type:
    sites <- which(site_type == type)
    spp_1 <- which(sp_group == type * 2 - 1)
    spp_2 <- which(sp_group == type * 2)
    
    # Create data.table with all interactions:
    data <- lapply(sites, function(x) {
      site_prop <- eval(parse(text = as.character(p_sample)))
      data <- expand.grid(
        site_id = x,
        sp1_id = sample(spp_1, max(round(length(spp_1) * site_prop), 1)), 
        sp2_id = sample(spp_2, max(round(length(spp_2) * site_prop), 1))
      ) %>% data.table::as.data.table()
    }) %>% data.table::rbindlist()
  }) %>% data.table::rbindlist()
  
  # Calculate standardized bioclimatic suitability product of interactions:
  data[, ':='(
    SS = scale(S[cbind(sp1_id, site_id)] * S[cbind(sp2_id, site_id)])
  )]
  
  # Compute response variable:
  data[, p := boot::inv.logit(
    alpha + beta[site_id] + gamma[sp1_id] + gamma[sp2_id] + lambda[site_id] * SS
  )]
  data[, Y := purrr::map_int(p, ~rbinom(1, 1, .))]
  
  # Remove empty rows/columns:
  if (rm_empty == TRUE) {
    data[, sum_sp1_ints := sum(Y), by = .(site_id, sp1_id)]
    data[, sum_sp2_ints := sum(Y), by = .(site_id, sp2_id)]
    data <- data[sum_sp1_ints > 0 & sum_sp2_ints > 0]
  }
  
  # Return data specified as a list:
  list(
    nb_sites = nb_sites,
    nb_spp = nb_spp,
    nb_int = nrow(data),
    nb_types = nb_types,
    nb_groups = 2 * nb_types,
    Y_array = data[, .(Y, site_id, sp1_id, sp2_id)],
    site_type = site_type,
    sp_group = sp_group,
    SS = data$SS,
    alpha = alpha,
    lambda_bar = lambda_bar,
    sigma_beta = sigma_beta,
    sigma_gamma = sigma_gamma,
    sigma_lambda = sigma_lambda,
    beta = beta,
    gamma = gamma,
    lambda = lambda
  )
}

generate_stan_start_values.all_binom_03 <- function(
  nb_sites, nb_spp, nb_types, p_sample =  "1.0", rm_empty = TRUE
) {
  list(
    alpha = 0,
    lambda_bar = rep(0, nb_types),
    zbeta = rep(0, nb_sites),
    zgamma = rep(0, nb_spp),
    zlambda = rep(0, nb_sites),
    sigma_beta = 0.1,
    sigma_gamma = rep(0.1, 2 * nb_types),
    sigma_lambda = rep(0.1, nb_types)
  )
}
