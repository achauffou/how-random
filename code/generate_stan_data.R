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


# Specific functions to generate data and starting values for pollination ======
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
  beta <- rnorm(nb_sites, 0, 1) %>% scale(scale = FALSE) %>% as.vector()
  gamma_pla <- rnorm(nb_pla, 0, 1) %>% scale(scale = FALSE) %>% as.vector()
  gamma_pol <- rnorm(nb_pol, 0, 1) %>% scale(scale = FALSE) %>% as.vector()
  lambda_bar <- rbeta(1, 2, 4)
  lambda <- rnorm(nb_sites, 0, 1) %>% scale(scale = FALSE) %>% 
    add(lambda_bar) %>% as.vector()
  
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
  
  # Update parameters and data if some species or sites were removed by filters:
  inc_sites <- sort(unique(data$site_id))
  inc_pla <- sort(unique(data$pla_id))
  inc_pol <- sort(unique(data$pol_id))
  beta <- beta[inc_sites]
  gamma_pla <- gamma_pla[inc_pla]
  gamma_pol <- gamma_pol[inc_pol]
  lambda <- lambda[inc_sites]
  data[, site_id := which(inc_sites == unique(site_id)), by = .(site_id)]
  data[, pla_id := which(inc_pla == unique(pla_id)), by = .(pla_id)]
  data[, pol_id := which(inc_pol == unique(pol_id)), by = .(pol_id)]
  
  # Return data specified as a list:
  list(
    nb_sites = length(inc_sites),
    nb_pla = length(inc_pla),
    nb_pol = length(inc_pol),
    nb_int = nrow(data),
    Y_array = data[, .(Y, site_id, pla_id, pol_id, rep = as.integer(1))],
    SS = data$SS,
    alpha = alpha + mean(beta) + mean(gamma_pla) + mean(gamma_pol),
    lambda_bar = mean(lambda_bar),
    sigma_beta = sd(beta),
    sigma_gamma_pla = sd(gamma_pla),
    sigma_gamma_pol = sd(gamma_pol),
    sigma_lambda = sd(lambda),
    beta = beta %>% scale(scale = FALSE) %>% as.vector(),
    gamma_pla = gamma_pla %>% scale(scale = FALSE) %>% as.vector(),
    gamma_pol = gamma_pol %>% scale(scale = FALSE) %>% as.vector(),
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

#' Generate data for pollination interactions with two suitability terms
#' 
generate_stan_data.pol_binom_09 <- function(
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
    S_pla = scale(S_pla),
    S_pol = scale(S_pol)
  )]
  
  # Sample parameters (center intercept parameters):
  alpha <- rbeta(1, 2, 4) * -1
  beta <- rnorm(nb_sites, 0, 1) %>% scale(scale = FALSE) %>% as.vector()
  gamma_pla <- rnorm(nb_pla, 0, 1) %>% scale(scale = FALSE) %>% as.vector()
  gamma_pol <- rnorm(nb_pol, 0, 1) %>% scale(scale = FALSE) %>% as.vector()
  lambda_pla <- rbeta(1, 2, 4)
  lambda_pol <- rbeta(1, 2, 4)
  
  # Compute response variable:
  data[, p := boot::inv.logit(
    alpha + beta[site_id] + gamma_pla[pla_id] + gamma_pol[pol_id] +
      lambda_pla * S_pla + lambda_pol * S_pol
  )]
  data[, Y := purrr::map_int(p, ~rbinom(1, 1, .))]
  
  # Remove empty rows/columns:
  if (rm_empty == TRUE) {
    data[, sum_pla_ints := sum(Y), by = .(site_id, pla_id)]
    data[, sum_pol_ints := sum(Y), by = .(site_id, pol_id)]
    data <- data[sum_pla_ints > 0 & sum_pol_ints > 0]
  }
  
  # Update parameters and data if some species or sites were removed by filters:
  inc_sites <- sort(unique(data$site_id))
  inc_pla <- sort(unique(data$pla_id))
  inc_pol <- sort(unique(data$pol_id))
  beta <- beta[inc_sites]
  gamma_pla <- gamma_pla[inc_pla]
  gamma_pol <- gamma_pol[inc_pol]
  data[, site_id := which(inc_sites == unique(site_id)), by = .(site_id)]
  data[, pla_id := which(inc_pla == unique(pla_id)), by = .(pla_id)]
  data[, pol_id := which(inc_pol == unique(pol_id)), by = .(pol_id)]
  
  # Return data specified as a list:
  list(
    nb_sites = length(inc_sites),
    nb_pla = length(inc_pla),
    nb_pol = length(inc_pol),
    nb_int = nrow(data),
    Y_array = data[, .(Y, site_id, pla_id, pol_id, rep = as.integer(1))],
    S_pla = data$S_pla,
    S_pol = data$S_pol,
    alpha = alpha + mean(beta) + mean(gamma_pla) + mean(gamma_pol),
    lambda_pla = lambda_pla,
    lambda_pol = lambda_pol,
    sigma_beta = sd(beta),
    sigma_gamma_pla = sd(gamma_pla),
    sigma_gamma_pol = sd(gamma_pol),
    beta = beta %>% scale(scale = FALSE) %>% as.vector(),
    gamma_pla = gamma_pla %>% scale(scale = FALSE) %>% as.vector(),
    gamma_pol = gamma_pol %>% scale(scale = FALSE) %>% as.vector()
  )
}

generate_stan_start_values.pol_binom_09 <- function(
  nb_sites, nb_pla, nb_pol, p_sample = "1.0", rm_empty = TRUE
) {
  list(
    alpha = 0,
    lambda_pla = 0,
    lambda_pol = 0,
    zbeta = rep(0, nb_sites),
    zgamma_pla = rep(0, nb_pla),
    zgamma_pol = rep(0, nb_pol),
    sigma_beta = 0.1,
    sigma_gamma_pla = 0.1,
    sigma_gamma_pol = 0.1
  )
}

#' Generate data for pollination binomial with origin status and single lambda
#' 
generate_stan_data.pol_binom_24 <- function(
  nb_sites, nb_pla, nb_pol, p_sample =  "1.0", rm_empty = TRUE
) {
  # Generate optimal suitability:
  S_opt_pla <- runif(nb_pla, 0, 1)
  S_opt_pol <- runif(nb_pol, 0, 1)
  env_sit <- runif(nb_sites, 0, 1)
  S_pla <- 1 - abs(outer(S_opt_pla, env_sit, "-"))
  S_pol <- 1 - abs(outer(S_opt_pol, env_sit, "-"))
  
  # Generate origin status:
  nat_pla <- matrix(rbinom(nb_pla * nb_sites, 1, 0.8), nrow = nb_pla, ncol = nb_sites)
  nat_pol <- matrix(rbinom(nb_pol * nb_sites, 1, 0.8), nrow = nb_pol, ncol = nb_sites)
  
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
    S_pol = S_pol[cbind(pol_id, site_id)],
    nat_pla = nat_pla[cbind(pla_id, site_id)],
    nat_pol = nat_pol[cbind(pol_id, site_id)]
  )]
  
  # Standardize predictor variables:
  data[, ':='(
    SS = scale(S_pla * S_pol)
  )]
  
  # Sample parameters (center intercept parameters):
  alpha <- rbeta(1, 2, 4) * -1
  beta <- rnorm(nb_sites, 0, 1) %>% scale(scale = FALSE) %>% as.vector()
  gamma_pla <- rnorm(nb_pla, 0, 1) %>% scale(scale = FALSE) %>% as.vector()
  gamma_pol <- rnorm(nb_pol, 0, 1) %>% scale(scale = FALSE) %>% as.vector()
  lambda <- rbeta(1, 2, 4)
  mu_pla <- rnorm(1, 0, 1)
  mu_pol <- rnorm(1, 0, 1)
  
  # Compute response variable:
  data[, p := boot::inv.logit(
    alpha + beta[site_id] + gamma_pla[pla_id] + gamma_pol[pol_id] +
      lambda * SS + nat_pla * mu_pla + nat_pol * mu_pol
  )]
  data[, Y := purrr::map_int(p, ~rbinom(1, 1, .))]
  
  # Remove empty rows/columns:
  if (rm_empty == TRUE) {
    data[, sum_pla_ints := sum(Y), by = .(site_id, pla_id)]
    data[, sum_pol_ints := sum(Y), by = .(site_id, pol_id)]
    data <- data[sum_pla_ints > 0 & sum_pol_ints > 0]
  }
  
  # Update parameters and data if some species or sites were removed by filters:
  inc_sites <- sort(unique(data$site_id))
  inc_pla <- sort(unique(data$pla_id))
  inc_pol <- sort(unique(data$pol_id))
  beta <- beta[inc_sites]
  gamma_pla <- gamma_pla[inc_pla]
  gamma_pol <- gamma_pol[inc_pol]
  data[, site_id := which(inc_sites == unique(site_id)), by = .(site_id)]
  data[, pla_id := which(inc_pla == unique(pla_id)), by = .(pla_id)]
  data[, pol_id := which(inc_pol == unique(pol_id)), by = .(pol_id)]
  
  # Return data specified as a list:
  list(
    nb_sites = length(inc_sites),
    nb_pla = length(inc_pla),
    nb_pol = length(inc_pol),
    nb_int = nrow(data),
    Y_array = data[, .(Y, site_id, pla_id, pol_id, rep = as.integer(1), nat_pla, nat_pol)],
    SS = data$SS,
    alpha = alpha + mean(beta) + mean(gamma_pla) + mean(gamma_pol),
    lambda = lambda,
    mu_pla = mu_pla,
    mu_pol = mu_pol,
    sigma_beta = sd(beta),
    sigma_gamma_pla = sd(gamma_pla),
    sigma_gamma_pol = sd(gamma_pol),
    beta = beta %>% scale(scale = FALSE) %>% as.vector(),
    gamma_pla = gamma_pla %>% scale(scale = FALSE) %>% as.vector(),
    gamma_pol = gamma_pol %>% scale(scale = FALSE) %>% as.vector(),
    SS_mean = mean(data$S_pla * data$S_pol),
    SS_sd = sd(data$S_pla * data$S_pol)
  )
}

generate_stan_start_values.pol_binom_24 <- function(
  nb_sites, nb_pla, nb_pol, p_sample = "1.0", rm_empty = TRUE
) {
  list(
    alpha = 0,
    lambda = 0,
    mu_pla = 0,
    mu_pol = 0,
    zbeta = rep(0, nb_sites),
    zgamma_pla = rep(0, nb_pla),
    zgamma_pol = rep(0, nb_pol),
    zlambda = rep(0, nb_sites),
    sigma_beta = 0.1,
    sigma_gamma_pla = 0.1,
    sigma_gamma_pol = 0.1
  )
}


# Specific functions to generate data and starting values for several ints =====
#' Generate data for all interaction types binomial with intercepts and slopes
#' 
generate_stan_data.all_binom_03 <- function(
  nb_sites, nb_spp, nb_types, p_sample =  "1.0", rm_empty = TRUE
) {
  # Assign all species and sites to an interaction type and group:
  site_type <- sample(1:nb_types, nb_sites, replace = TRUE) %>% sort()
  sp_group <- sample(1:(2 * nb_types), nb_spp, replace = TRUE) %>% sort()
  
  # Sample parameters:
  alpha <- rbeta(nb_types, 4, 2) * -1.5 - 1.0
  lambda_bar <- rbeta(nb_types, 4, 2) + 0.1
  sigma_beta <- rbeta(nb_types, 2, 3) + 0.5
  sigma_gamma <- rbeta(2 * nb_types, 2, 3) + 0.5
  sigma_lambda <- rbeta(nb_types, 2, 3) + 0.5
  zbeta <- rnorm(nb_sites, 0, 1) %>% scale() %>% as.vector()
  zgamma <- rnorm(nb_spp, 0, 1) %>% scale() %>% as.vector()
  zlambda <- rnorm(nb_sites, 0, 1) %>% scale() %>% as.vector()
  
  # Compute non-centered parametrization:
  beta <- zbeta * sigma_beta[site_type]
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
      if (length(spp_1) > 1) {
        the_spp_1 <- sample(spp_1, max(round(length(spp_1) * site_prop), 1))
      } else {
        the_spp_1 <- spp_1
      }
      if (length(spp_2) > 1) {
        the_spp_2 <- sample(spp_2, max(round(length(spp_2) * site_prop), 1))
      } else {
        the_spp_2 <- spp_2
      }
      data <- expand.grid(
        site_id = x,
        sp1_id = the_spp_1, 
        sp2_id = the_spp_2
      ) %>% data.table::as.data.table()
    }) %>% data.table::rbindlist()
  }) %>% data.table::rbindlist()
  
  # Calculate standardized bioclimatic suitability product of interactions:
  data[, ':='(
    SS = scale(S[cbind(sp1_id, site_id)] * S[cbind(sp2_id, site_id)])
  )]
  
  # Compute response variable:
  data[, p := boot::inv.logit(
    alpha[site_type[site_id]] + beta[site_id] + gamma[sp1_id] + gamma[sp2_id] + 
      lambda[site_id] * SS
  )]
  data[, Y := purrr::map_int(p, ~rbinom(1, 1, .))]
  
  # Remove empty rows/columns:
  if (rm_empty == TRUE) {
    data[, sum_sp1_ints := sum(Y), by = .(site_id, sp1_id)]
    data[, sum_sp2_ints := sum(Y), by = .(site_id, sp2_id)]
    data <- data[sum_sp1_ints > 0 & sum_sp2_ints > 0]
  }
  
  # Update parameters and data if some species or sites were removed by filters:
  inc_sites <- sort(unique(data$site_id))
  inc_spp <- sort(unique(c(data$sp1_id, data$sp2_id)))
  data[, site_id := which(inc_sites == unique(site_id)), by = .(site_id)]
  data[, sp1_id := which(inc_spp == unique(sp1_id)), by = .(sp1_id)]
  data[, sp2_id := which(inc_spp == unique(sp2_id)), by = .(sp2_id)]
  zbeta <- zbeta[inc_sites]
  zgamma <- zgamma[inc_spp]
  beta <- beta[inc_sites]
  gamma <- gamma[inc_spp]
  lambda <- lambda[inc_sites]
  site_type <- site_type[inc_sites]
  sp_group <- sp_group[inc_spp]
  for (int in 1:nb_types) {
    sigma_beta[int] <- sd(beta[site_type == int])
    alpha[int] <- alpha[int] + mean(beta[site_type == int]) + 
      mean(gamma[sp_group == int * 2 - 1]) + mean(gamma[sp_group == int * 2])
    if (length(gamma[sp_group == 2 * int - 1]) > 1) {
      sigma_gamma[2 * int - 1] <- sd(gamma[sp_group == 2 * int - 1])
    } else {
      sigma_gamma[2 * int - 1] <- -1
    }
    if (length(gamma[sp_group == 2 * int]) > 1) {
      sigma_gamma[2 * int] <- sd(gamma[sp_group == 2 * int])
    } else {
      sigma_gamma[2 * int] <- -1
    }
    sigma_lambda[int] <- sd(lambda[site_type == int])
    lambda_bar[int] <- mean(lambda[site_type == int])
  }
  zbeta <- zbeta %>% scale() %>% as.vector()
  zgamma <-zgamma %>% scale() %>% as.vector()
  zlambda <- zbeta %>% scale() %>% as.vector()
  beta <- zbeta * sigma_beta
  gamma <- zgamma * sigma_gamma[sp_group]
  lambda <- zlambda * sigma_lambda[site_type] + lambda_bar[site_type]
  
  # Return data specified as a list:
  list(
    nb_sites = length(inc_sites),
    nb_spp = length(inc_spp),
    nb_int = nrow(data),
    nb_types = nb_types,
    nb_groups = 2 * nb_types,
    Y_array = data[, .(Y, site_id, sp1_id, sp2_id, rep = as.integer(1))],
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
    alpha = rep(0, nb_types),
    lambda_bar = rep(0, nb_types),
    zbeta = rep(0, nb_sites),
    zgamma = rep(0, nb_spp),
    zlambda = rep(0, nb_sites),
    sigma_beta = rep(0.1, nb_types),
    sigma_gamma = rep(0.1, 2 * nb_types),
    sigma_lambda = rep(0.1, nb_types)
  )
}

#' Generate data for all interaction types binomial with intercepts and slope
#' 
generate_stan_data.all_binom_04 <- function(
  nb_sites, nb_spp, nb_types, p_sample =  "1.0", rm_empty = TRUE
) {
  # Assign all species and sites to an interaction type and group:
  site_type <- sample(1:nb_types, nb_sites, replace = TRUE) %>% sort()
  sp_group <- sample(1:(2 * nb_types), nb_spp, replace = TRUE) %>% sort()
  
  # Sample parameters:
  alpha <- rbeta(nb_types, 4, 2) * -1.5 - 1.0
  lambda <- rbeta(nb_types, 4, 2) - 0.1
  sigma_beta <- rbeta(nb_types, 2, 3) + 0.5
  sigma_gamma <- rbeta(2 * nb_types, 2, 3) + 0.5
  zbeta <- rnorm(nb_sites, 0, 1) %>% scale() %>% as.vector()
  zgamma <- rnorm(nb_spp, 0, 1) %>% scale() %>% as.vector()
  
  # Compute non-centered parametrization:
  beta <- zbeta * sigma_beta[site_type]
  gamma <- zgamma * sigma_gamma[sp_group]
  
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
      if (length(spp_1) > 1) {
        the_spp_1 <- sample(spp_1, max(round(length(spp_1) * site_prop), 1))
      } else {
        the_spp_1 <- spp_1
      }
      if (length(spp_2) > 1) {
        the_spp_2 <- sample(spp_2, max(round(length(spp_2) * site_prop), 1))
      } else {
        the_spp_2 <- spp_2
      }
      data <- expand.grid(
        site_id = x,
        sp1_id = the_spp_1, 
        sp2_id = the_spp_2
      ) %>% data.table::as.data.table()
    }) %>% data.table::rbindlist()
  }) %>% data.table::rbindlist()
  
  # Calculate standardized bioclimatic suitability product of interactions:
  data[, ':='(
    SS = scale(S[cbind(sp1_id, site_id)] * S[cbind(sp2_id, site_id)])
  )]
  
  # Compute response variable:
  data[, p := boot::inv.logit(
    alpha[site_type[site_id]] + beta[site_id] + gamma[sp1_id] + gamma[sp2_id] + 
      lambda[site_type[site_id]] * SS
  )]
  data[, Y := purrr::map_int(p, ~rbinom(1, 1, .))]
  
  # Remove empty rows/columns:
  if (rm_empty == TRUE) {
    data[, sum_sp1_ints := sum(Y), by = .(site_id, sp1_id)]
    data[, sum_sp2_ints := sum(Y), by = .(site_id, sp2_id)]
    data <- data[sum_sp1_ints > 0 & sum_sp2_ints > 0]
  }
  
  # Update parameters and data if some species or sites were removed by filters:
  inc_sites <- sort(unique(data$site_id))
  inc_spp <- sort(unique(c(data$sp1_id, data$sp2_id)))
  data[, site_id := which(inc_sites == unique(site_id)), by = .(site_id)]
  data[, sp1_id := which(inc_spp == unique(sp1_id)), by = .(sp1_id)]
  data[, sp2_id := which(inc_spp == unique(sp2_id)), by = .(sp2_id)]
  zbeta <- zbeta[inc_sites]
  zgamma <- zgamma[inc_spp]
  beta <- beta[inc_sites]
  gamma <- gamma[inc_spp]
  site_type <- site_type[inc_sites]
  sp_group <- sp_group[inc_spp]
  for (int in 1:nb_types) {
    sigma_beta[int] <- sd(beta[site_type == int])
    alpha[int] <- alpha[int] + mean(beta[site_type == int]) + 
      mean(gamma[sp_group == int * 2 - 1]) + mean(gamma[sp_group == int * 2])
    if (length(gamma[sp_group == 2 * int - 1]) > 1) {
      sigma_gamma[2 * int - 1] <- sd(gamma[sp_group == 2 * int - 1])
    } else {
      sigma_gamma[2 * int - 1] <- -1
    }
    if (length(gamma[sp_group == 2 * int]) > 1) {
      sigma_gamma[2 * int] <- sd(gamma[sp_group == 2 * int])
    } else {
      sigma_gamma[2 * int] <- -1
    }
  }
  zbeta <- zbeta %>% scale() %>% as.vector()
  zgamma <-zgamma %>% scale() %>% as.vector()
  beta <- zbeta * sigma_beta
  gamma <- zgamma * sigma_gamma[sp_group]
  
  # Return data specified as a list:
  list(
    nb_sites = length(inc_sites),
    nb_spp = length(inc_spp),
    nb_int = nrow(data),
    nb_types = nb_types,
    nb_groups = 2 * nb_types,
    Y_array = data[, .(Y, site_id, sp1_id, sp2_id, rep = as.integer(1))],
    site_type = site_type,
    sp_group = sp_group,
    SS = data$SS,
    alpha = alpha,
    lambda = lambda,
    sigma_beta = sigma_beta,
    sigma_gamma = sigma_gamma,
    beta = beta,
    gamma = gamma
  )
}

generate_stan_start_values.all_binom_04 <- function(
  nb_sites, nb_spp, nb_types, p_sample =  "1.0", rm_empty = TRUE
) {
  list(
    alpha = rep(0, nb_types),
    lambda_bar = rep(0, nb_types),
    zbeta = rep(0, nb_sites),
    zgamma = rep(0, nb_spp),
    sigma_beta = rep(0.1, nb_types),
    sigma_gamma = rep(0.1, 2 * nb_types)
  )
}

#' Generate data for all interaction types binomial with two suitability terms
#' 
generate_stan_data.all_binom_09 <- function(
  nb_sites, nb_spp, nb_types, p_sample =  "1.0", rm_empty = TRUE
) {
  # Assign all species and sites to an interaction type and group:
  site_type <- sample(1:nb_types, nb_sites, replace = TRUE) %>% sort()
  sp_group <- sample(1:(2 * nb_types), nb_spp, replace = TRUE) %>% sort()
  
  # Sample parameters:
  alpha <- rbeta(nb_types, 4, 2) * -1.5 - 1.0
  lambda <- rbeta(nb_types * 2, 4, 2) - 0.1
  sigma_beta <- rbeta(nb_types, 2, 3) + 0.5
  sigma_gamma <- rbeta(2 * nb_types, 2, 3) + 0.5
  zbeta <- rnorm(nb_sites, 0, 1) %>% scale() %>% as.vector()
  zgamma <- rnorm(nb_spp, 0, 1) %>% scale() %>% as.vector()
  
  # Compute non-centered parametrization:
  beta <- zbeta * sigma_beta[site_type]
  gamma <- zgamma * sigma_gamma[sp_group]
  
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
      if (length(spp_1) > 1) {
        the_spp_1 <- sample(spp_1, max(round(length(spp_1) * site_prop), 1))
      } else {
        the_spp_1 <- spp_1
      }
      if (length(spp_2) > 1) {
        the_spp_2 <- sample(spp_2, max(round(length(spp_2) * site_prop), 1))
      } else {
        the_spp_2 <- spp_2
      }
      data <- expand.grid(
        site_id = x,
        sp1_id = the_spp_1, 
        sp2_id = the_spp_2
      ) %>% data.table::as.data.table()
    }) %>% data.table::rbindlist()
  }) %>% data.table::rbindlist()
  
  # Calculate standardized bioclimatic suitability product of interactions:
  data[, ':='(
    sp1_group = sp_group[sp1_id],
    sp2_group = sp_group[sp2_id],
    S1 = S[cbind(sp1_id, site_id)],
    S2 = S[cbind(sp2_id, site_id)]
  )]
  S1_mean <- data[, .(S1_mean = mean(S1)), by = .(sp1_group)][order(sp1_group)][['S1_mean']]
  S1_sd <- data[, .(S1_sd = sd(S1)), by = .(sp1_group)][order(sp1_group)][['S1_sd']]
  S2_mean <- data[, .(S2_mean = mean(S2)), by = .(sp2_group)][order(sp2_group)][['S2_mean']]
  S2_sd <- data[, .(S2_sd = sd(S2)), by = .(sp2_group)][order(sp2_group)][['S2_sd']]
  S_mean <- numeric(2 * nb_types)
  S_sd <- numeric(2 * nb_types)
  for (i in 1:nb_types) {
    S_mean[2 * i - 1] <- S1_mean[i]
    S_mean[2 * i] <- S2_mean[i]
    S_sd[2 * i - 1] <- S1_sd[i]
    S_sd[2 * i] <- S2_sd[i]
  }
  data[, ':='(
    S1 = (S1 - S_mean[sp1_group]) / S_sd[sp1_group],
    S2 = (S2 - S_mean[sp2_group]) / S_sd[sp2_group]
  )]
  
  # Compute response variable:
  data[, p := boot::inv.logit(
    alpha[site_type[site_id]] + beta[site_id] + gamma[sp1_id] + gamma[sp2_id] + 
      lambda[sp_group[sp1_id]] * S1 + lambda[sp_group[sp2_id]] * S2
  )]
  data[, Y := purrr::map_int(p, ~rbinom(1, 1, .))]
  
  # Remove empty rows/columns:
  if (rm_empty == TRUE) {
    data[, sum_sp1_ints := sum(Y), by = .(site_id, sp1_id)]
    data[, sum_sp2_ints := sum(Y), by = .(site_id, sp2_id)]
    data <- data[sum_sp1_ints > 0 & sum_sp2_ints > 0]
  }
  
  # Update parameters and data if some species or sites were removed by filters:
  inc_sites <- sort(unique(data$site_id))
  inc_spp <- sort(unique(c(data$sp1_id, data$sp2_id)))
  data[, site_id := which(inc_sites == unique(site_id)), by = .(site_id)]
  data[, sp1_id := which(inc_spp == unique(sp1_id)), by = .(sp1_id)]
  data[, sp2_id := which(inc_spp == unique(sp2_id)), by = .(sp2_id)]
  zbeta <- zbeta[inc_sites]
  zgamma <- zgamma[inc_spp]
  beta <- beta[inc_sites]
  gamma <- gamma[inc_spp]
  site_type <- site_type[inc_sites]
  sp_group <- sp_group[inc_spp]
  for (int in 1:nb_types) {
    sigma_beta[int] <- sd(beta[site_type == int])
    alpha[int] <- alpha[int] + mean(beta[site_type == int]) + 
      mean(gamma[sp_group == int * 2 - 1]) + mean(gamma[sp_group == int * 2])
    if (length(gamma[sp_group == 2 * int - 1]) > 1) {
      sigma_gamma[2 * int - 1] <- sd(gamma[sp_group == 2 * int - 1])
    } else {
      sigma_gamma[2 * int - 1] <- -1
    }
    if (length(gamma[sp_group == 2 * int]) > 1) {
      sigma_gamma[2 * int] <- sd(gamma[sp_group == 2 * int])
    } else {
      sigma_gamma[2 * int] <- -1
    }
  }
  zbeta <- zbeta %>% scale() %>% as.vector()
  zgamma <-zgamma %>% scale() %>% as.vector()
  beta <- zbeta * sigma_beta
  gamma <- zgamma * sigma_gamma[sp_group]
  
  # Return data specified as a list:
  list(
    nb_sites = length(inc_sites),
    nb_spp = length(inc_spp),
    nb_int = nrow(data),
    nb_types = nb_types,
    nb_groups = 2 * nb_types,
    Y_array = data[, .(Y, site_id, sp1_id, sp2_id, rep = as.integer(1))],
    site_type = site_type,
    sp_group = sp_group,
    S1 = data$S1,
    S2 = data$S2,
    S_mean = S_mean,
    S_sd = S_sd,
    alpha = alpha,
    lambda = lambda,
    sigma_beta = sigma_beta,
    sigma_gamma = sigma_gamma,
    beta = beta,
    gamma = gamma
  )
}

generate_stan_start_values.all_binom_09 <- function(
  nb_sites, nb_spp, nb_types, p_sample =  "1.0", rm_empty = TRUE
) {
  list(
    alpha = rep(0, nb_types),
    lambda = rep(0, 2 * nb_types),
    zbeta = rep(0, nb_sites),
    zgamma = rep(0, nb_spp),
    sigma_beta = rep(0.1, nb_types),
    sigma_gamma = rep(0.1, 2 * nb_types)
  )
}

#' Generate data for all interaction types binomial with origin status
#' 
generate_stan_data.all_binom_24 <- function(
  nb_sites, nb_spp, nb_types, p_sample =  "1.0", rm_empty = TRUE
) {
  # Assign all species and sites to an interaction type and group:
  site_type <- sample(1:nb_types, nb_sites, replace = TRUE) %>% sort()
  sp_group <- sample(1:(2 * nb_types), nb_spp, replace = TRUE) %>% sort()
  
  # Sample parameters:
  alpha <- rbeta(nb_types, 4, 2) * -1.5 - 1.0
  lambda <- rbeta(nb_types, 4, 2) - 0.1
  mu <- - rbeta(2 * nb_types, 2, 3) - 0.25
  sigma_beta <- rbeta(nb_types, 2, 3) + 0.5
  sigma_gamma <- rbeta(2 * nb_types, 2, 3) + 0.5
  zbeta <- rnorm(nb_sites, 0, 1) %>% scale() %>% as.vector()
  zgamma <- rnorm(nb_spp, 0, 1) %>% scale() %>% as.vector()
  
  # Compute non-centered parametrization:
  beta <- zbeta * sigma_beta[site_type]
  gamma <- zgamma * sigma_gamma[sp_group]
  
  # Generate optimal suitability:
  S_opt <- runif(nb_spp, 0, 1)
  env_sit <- runif(nb_sites, 0, 1)
  S <- 1 - abs(outer(S_opt, env_sit, "-"))
  
  # Generate origin status:
  nat <- matrix(rbinom(nb_spp * nb_sites, 1, 0.8), nrow = nb_spp, ncol = nb_sites)
  
  # Generate data interaction type by interaction type:
  data <- lapply(1:nb_types, function(type) {
    # Get IDs of species and sites that belong to this interaction type:
    sites <- which(site_type == type)
    spp_1 <- which(sp_group == type * 2 - 1)
    spp_2 <- which(sp_group == type * 2)
    
    # Create data.table with all interactions:
    data <- lapply(sites, function(x) {
      site_prop <- eval(parse(text = as.character(p_sample)))
      if (length(spp_1) > 1) {
        the_spp_1 <- sample(spp_1, max(round(length(spp_1) * site_prop), 1))
      } else {
        the_spp_1 <- spp_1
      }
      if (length(spp_2) > 1) {
        the_spp_2 <- sample(spp_2, max(round(length(spp_2) * site_prop), 1))
      } else {
        the_spp_2 <- spp_2
      }
      data <- expand.grid(
        site_id = x,
        sp1_id = the_spp_1, 
        sp2_id = the_spp_2
      ) %>% data.table::as.data.table()
    }) %>% data.table::rbindlist()
  }) %>% data.table::rbindlist()
  
  # Calculate standardized bioclimatic suitability product of interactions:
  data[, ':='(
    SS = scale(S[cbind(sp1_id, site_id)] * S[cbind(sp2_id, site_id)]),
    sp1_is_native = nat[cbind(sp1_id, site_id)],
    sp2_is_native = nat[cbind(sp2_id, site_id)]
  )]
  
  # Compute response variable:
  data[, p := boot::inv.logit(
    alpha[site_type[site_id]] + beta[site_id] + gamma[sp1_id] + gamma[sp2_id] + 
      lambda[site_type[site_id]] * SS + sp1_is_native * mu[sp_group[sp1_id]] +
      sp2_is_native * mu[sp_group[sp2_id]]
  )]
  data[, Y := purrr::map_int(p, ~rbinom(1, 1, .))]
  
  # Remove empty rows/columns:
  if (rm_empty == TRUE) {
    data[, sum_sp1_ints := sum(Y), by = .(site_id, sp1_id)]
    data[, sum_sp2_ints := sum(Y), by = .(site_id, sp2_id)]
    data <- data[sum_sp1_ints > 0 & sum_sp2_ints > 0]
  }
  
  # Update parameters and data if some species or sites were removed by filters:
  inc_sites <- sort(unique(data$site_id))
  inc_spp <- sort(unique(c(data$sp1_id, data$sp2_id)))
  data[, site_id := which(inc_sites == unique(site_id)), by = .(site_id)]
  data[, sp1_id := which(inc_spp == unique(sp1_id)), by = .(sp1_id)]
  data[, sp2_id := which(inc_spp == unique(sp2_id)), by = .(sp2_id)]
  zbeta <- zbeta[inc_sites]
  zgamma <- zgamma[inc_spp]
  beta <- beta[inc_sites]
  gamma <- gamma[inc_spp]
  site_type <- site_type[inc_sites]
  sp_group <- sp_group[inc_spp]
  for (int in 1:nb_types) {
    sigma_beta[int] <- sd(beta[site_type == int])
    alpha[int] <- alpha[int] + mean(beta[site_type == int]) + 
      mean(gamma[sp_group == int * 2 - 1]) + mean(gamma[sp_group == int * 2])
    if (length(gamma[sp_group == 2 * int - 1]) > 1) {
      sigma_gamma[2 * int - 1] <- sd(gamma[sp_group == 2 * int - 1])
    } else {
      sigma_gamma[2 * int - 1] <- -1
    }
    if (length(gamma[sp_group == 2 * int]) > 1) {
      sigma_gamma[2 * int] <- sd(gamma[sp_group == 2 * int])
    } else {
      sigma_gamma[2 * int] <- -1
    }
  }
  zbeta <- zbeta %>% scale() %>% as.vector()
  zgamma <-zgamma %>% scale() %>% as.vector()
  beta <- zbeta * sigma_beta
  gamma <- zgamma * sigma_gamma[sp_group]
  
  # Return data specified as a list:
  list(
    nb_sites = length(inc_sites),
    nb_spp = length(inc_spp),
    nb_int = nrow(data),
    nb_types = nb_types,
    nb_groups = 2 * nb_types,
    Y_array = data[, .(Y, site_id, sp1_id, sp2_id, rep = as.integer(1), sp1_is_native, sp2_is_native)],
    site_type = site_type,
    sp_group = sp_group,
    SS = data$SS,
    alpha = alpha,
    lambda = lambda,
    mu = mu,
    sigma_beta = sigma_beta,
    sigma_gamma = sigma_gamma,
    beta = beta,
    gamma = gamma
  )
}

generate_stan_start_values.all_binom_24 <- function(
  nb_sites, nb_spp, nb_types, p_sample =  "1.0", rm_empty = TRUE
) {
  list(
    alpha = rep(0, nb_types),
    lambda_bar = rep(0, nb_types),
    mu = rep(0, nb_types * 2),
    zbeta = rep(0, nb_sites),
    zgamma = rep(0, nb_spp),
    sigma_beta = rep(0.1, nb_types),
    sigma_gamma = rep(0.1, 2 * nb_types)
  )
}
