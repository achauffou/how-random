# Overall functions to prepare Stan data and starting values ===================
#' Prepare data for a Stan model from its specification
#' 
prepare_stan_data <- function(spec, interactions, results_folder = "results/analyses") {
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
    fun_args <- list(interactions)
  } else {
    fun_args <- c(list(interactions), spec$data_generation_args)
  }
  paste0("prepare_stan_data.", spec$data_generation_model) %>%
    get() %>%
    do.call(args = fun_args) %>%
    saveRDS(file = out_path)
  out_path
}

#' Generate starting values for a Stan model from its specification
#' 
prepare_stan_start_values <- function(spec, interactions, results_folder = "results/analyses") {
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
    fun_args <- list(interactions)
  } else {
    fun_args <- c(list(interactions), spec$data_generation_args)
  }
  paste0("prepare_stan_start_values.", spec$data_generation_model) %>%
    get() %>%
    do.call(args = fun_args) %>%
    saveRDS(file = out_path)
  out_path
}


# Functions to prepare data for specific Stan analyses =========================
#' Pollination binomial model with only intercepts (no lambdas)
#' 
prepare_stan_data.pol_binom_02 <- function(interactions) {
  # Select relevant columns with complete cases:
  ints <- interactions[int_type == "Pollination"]
  ints <- ints[, .(
    int_strength, net_id, sp1_id, sp1_name, sp1_kingdom, sp2_id, sp2_name, 
    sp2_kingdom
  )]
  ints <- ints[complete.cases(ints[, -c("sp1_kingdom", "sp2_kingdom")])]
  
  # Use binary interactions:
  ints[, ':='(Y = ifelse(int_strength > 0, 1, 0))]
  
  # Merge replicates together:
  ints[, ':='(Y = sum(Y), N = .N), by = .(net_id, sp1_id, sp2_id)]
  ints <- unique(ints, by = c("net_id", "sp1_id", "sp2_id"))
  
  # Prepare sites, plant and pollinators IDs:
  ints[, pla_id := .GRP, by = .(sp1_id)]
  ints[, pol_id := .GRP, by = .(sp2_id)]
  ints[, site_id := .GRP, by = .(net_id)]
  pla_ids <- ints[, .(out = unique(sp1_id)), by = pla_id][order(pla_id)][['out']]
  pla_names <- ints[, .(out = unique(sp1_name)), by = pla_id][order(pla_id)][['out']]
  pol_ids <- ints[, .(out = unique(sp2_id)), by = pol_id][order(pol_id)][['out']]
  pol_names <- ints[, .(out = unique(sp2_name)), by = pol_id][order(pol_id)][['out']]
  site_names <- ints[, .(out = unique(net_id)), by = site_id][order(site_id)][['out']]
  
  # Return list of data:
  list(
    nb_sites = length(unique(ints$site_id)),
    nb_pla = length(unique(ints$pla_id)),
    nb_pol = length(unique(ints$pol_id)),
    nb_int = nrow(ints),
    Y_array = ints[, .(Y, site_id, pla_id, pol_id, N)],
    site_names = site_names,
    pla_ids = pla_ids,
    pla_names = pla_names,
    pol_ids = pol_ids,
    pol_names = pol_names
  )
}

prepare_stan_start_values.pol_binom_02 <- function(interactions) {
  data <- prepare_stan_data.pol_binom_02(interactions)
  list(
    alpha = 0,
    lambda_bar = 0,
    zbeta = rep(0, data$nb_sites),
    zgamma_pla = rep(0, data$nb_pla),
    zgamma_pol = rep(0, data$nb_pol),
    sigma_beta = 0.1,
    sigma_gamma_pla = 0.1,
    sigma_gamma_pol = 0.1
  )
}

#' Pollination binomial model with bioclimatic variables
#' 
prepare_stan_data.pol_binom_bioclim <- function(
  interactions, min_bioclim_occs = NULL, collec = FALSE
) {
  # Select relevant columns with complete cases:
  ints <- interactions[int_type == "Pollination"]
  if (!is.null(min_bioclim_occs)) {
    ints <- ints[sp1_nb_bioclim >= min_bioclim_occs & 
                   sp2_nb_bioclim >= min_bioclim_occs]
  }
  if (collec == TRUE) {
    ints <- ints[, .(
      int_strength, net_id, sp1_id, sp1_name, sp1_kingdom, 
      sp1_suitability = sp1_collec_suitability, sp2_id, sp2_name, sp2_kingdom, 
      sp2_suitability = sp2_collec_suitability
    )]
  } else {
    ints <- ints[, .(
      int_strength, net_id, sp1_id, sp1_name, sp1_kingdom, 
      sp1_suitability = sp1_indiv_suitability, sp2_id, sp2_name, sp2_kingdom, 
      sp2_suitability = sp2_indiv_suitability
    )]
  }
  ints <- ints[complete.cases(ints[, -c("sp1_kingdom", "sp2_kingdom")])]
  
  # Use binary interactions:
  ints[, ':='(Y = ifelse(int_strength > 0, 1, 0))]
  
  # Merge replicates together:
  ints[, ':='(Y = sum(Y), N = .N), by = .(net_id, sp1_id, sp2_id)]
  ints <- unique(ints, by = c("net_id", "sp1_id", "sp2_id"))
  
  # Remove sites and species that do not have both zeros and ones:
  ints[, ones_sites := sum(Y), by = .(net_id)]
  ints[, zeros_sites := sum(Y == 0), by = .(net_id)]
  ints[, ones_sp1 := sum(Y), by = .(sp1_id)]
  ints[, zeros_sp1 := sum(Y == 0), by = .(sp1_id)]
  ints[, ones_sp2 := sum(Y), by = .(sp2_id)]
  ints[, zeros_sp2 := sum(Y == 0), by = .(sp2_id)]
  ints <- ints[ones_sites > 0 & zeros_sites > 0 & ones_sp1 > 0 & 
                 zeros_sp1 > 0 & ones_sp2 > 0 & zeros_sp2 > 0]
  
  # Prepare sites, plant and pollinators IDs:
  ints[, pla_id := .GRP, by = .(sp1_id)]
  ints[, pol_id := .GRP, by = .(sp2_id)]
  ints[, site_id := .GRP, by = .(net_id)]
  pla_ids <- ints[, .(out = unique(sp1_id)), by = pla_id][order(pla_id)][['out']]
  pla_names <- ints[, .(out = unique(sp1_name)), by = pla_id][order(pla_id)][['out']]
  pol_ids <- ints[, .(out = unique(sp2_id)), by = pol_id][order(pol_id)][['out']]
  pol_names <- ints[, .(out = unique(sp2_name)), by = pol_id][order(pol_id)][['out']]
  site_names <- ints[, .(out = unique(net_id)), by = site_id][order(site_id)][['out']]
  
  # Compute standardized product of suitabilities:
  SS_mean <- mean(ints[['sp1_suitability']] * ints[['sp2_suitability']])
  SS_sd <- sd(ints[['sp1_suitability']] * ints[['sp2_suitability']])
  ints[, SS := scale(sp1_suitability * sp2_suitability)]
  
  # Return list of data:
  list(
    nb_sites = length(unique(ints$site_id)),
    nb_pla = length(unique(ints$pla_id)),
    nb_pol = length(unique(ints$pol_id)),
    nb_int = nrow(ints),
    Y_array = ints[, .(Y, site_id, pla_id, pol_id, N)],
    SS = ints$SS,
    SS_mean = SS_mean,
    SS_sd = SS_sd,
    site_names = site_names,
    pla_ids = pla_ids,
    pla_names = pla_names,
    pol_ids = pol_ids,
    pol_names = pol_names
  )
}

prepare_stan_start_values.pol_binom_bioclim <- function(
  interactions, min_bioclim_occs = NULL, collec = FALSE
) {
  data <- prepare_stan_data.pol_binom_bioclim(interactions, min_bioclim_occs, collec)
  list(
    alpha = 0,
    lambda_bar = 0,
    zbeta = rep(0, data$nb_sites),
    zgamma_pla = rep(0, data$nb_pla),
    zgamma_pol = rep(0, data$nb_pol),
    zlambda = rep(0, data$nb_sites),
    sigma_beta = 0.1,
    sigma_gamma_pla = 0.1,
    sigma_gamma_pol = 0.1,
    sigma_lambda = 0.1
  )
}

#' All interactions binomial model with bioclimatic variables
#' 
prepare_stan_data.all_binom_bioclim <- function(
  interactions, min_nb_ints = 100, min_bioclim_occs = NULL, collec = FALSE
) {
  ints <- interactions
  # Select relevant columns with complete cases:
  if (!is.null(min_bioclim_occs)) {
    ints <- ints[sp1_nb_bioclim >= min_bioclim_occs & 
                   sp2_nb_bioclim >= min_bioclim_occs]
  }
  if (collec == TRUE) {
    ints <- ints[, .(
      int_strength, int_type, net_id, sp1_fun_group, sp1ID = sp1_id, sp1_name, 
      sp1_kingdom, sp1_suitability = sp1_collec_suitability, sp2_fun_group, 
      sp2ID = sp2_id, sp2_name, sp2_kingdom, sp2_suitability = sp2_collec_suitability
    )]
  } else {
    ints <- ints[, .(
      int_strength, int_type, net_id, sp1_fun_group, sp1ID = sp1_id, sp1_name, 
      sp1_kingdom, sp1_suitability = sp1_indiv_suitability, sp2_fun_group, 
      sp2ID = sp2_id, sp2_name, sp2_kingdom, sp2_suitability = sp2_indiv_suitability
    )]
  }
  ints <- ints[complete.cases(ints[, -c("sp1_kingdom", "sp2_kingdom")])]
  
  # Use binary interactions:
  ints[, ':='(Y = ifelse(int_strength > 0, 1, 0))]
  
  # Merge replicates together:
  ints[, ':='(Y = sum(Y), N = .N), by = .(net_id, sp1ID, sp2ID)]
  ints <- unique(ints, by = c("net_id", "sp1ID", "sp2ID"))
  
  # Remove sites and species that do not have both zeros and ones:
  ints[, ones_sites := sum(Y), by = .(net_id)]
  ints[, zeros_sites := sum(Y == 0), by = .(net_id)]
  ints[, ones_sp1 := sum(Y), by = .(sp1ID)]
  ints[, zeros_sp1 := sum(Y == 0), by = .(sp1ID)]
  ints[, ones_sp2 := sum(Y), by = .(sp2ID)]
  ints[, zeros_sp2 := sum(Y == 0), by = .(sp2ID)]
  ints <- ints[ones_sites > 0 & zeros_sites > 0 & ones_sp1 > 0 & 
                 zeros_sp1 > 0 & ones_sp2 > 0 & zeros_sp2 > 0]
  
  # Use only interaction types with enough data:
  type_names <- ints[, .N, by = .(int_type)][N > min_nb_ints][["int_type"]]
  ints <- ints[int_type %in% type_names]
  
  # Prepare species IDs:
  ints[, sp1_id := .GRP, by = .(sp1ID, int_type, sp1_fun_group)]
  sp1_ids <- ints[, .(out = unique(sp1ID)), by = sp1_id][order(sp1_id)][['out']]
  sp1_names <- ints[, .(out = unique(sp1_name)), by = sp1_id][order(sp1_id)][['out']]
  ints[, sp2_id := .GRP + length(sp1_ids), by = .(sp2ID, int_type, sp2_fun_group)]
  sp2_ids <- ints[, .(out = unique(sp2ID)), by = sp2_id][order(sp2_id)][['out']]
  sp2_names <- ints[, .(out = unique(sp2_name)), by = sp2_id][order(sp2_id)][['out']]
  sp_ids <- c(sp1_ids, sp2_ids)
  sp_names <- c(sp1_names, sp2_names)
  
  # Prepare sites IDs:
  ints[, site_id := .GRP, by = .(net_id)]
  site_names <- ints[, .(out = unique(net_id)), by = site_id][order(site_id)][['out']]
  
  # Prepare interaction types and functional groups IDs:
  ints[, type_id := .GRP, by = .(int_type)]
  type_names <- ints[, .(out = unique(int_type)), by = type_id][order(type_id)][['out']]
  ints[, group1_id := type_id]
  ints[, group2_id := length(type_names) + type_id]
  group1_names <- ints[, .(out = unique(sp1_fun_group)), by = group1_id][order(group1_id)][['out']]
  group2_names <- ints[, .(out = unique(sp2_fun_group)), by = group2_id][order(group2_id)][['out']]
  group_names <- c(group1_names, group2_names)
  group_types <- c(type_names, type_names)
  
  # Get the interaction type of each site and functional group of each species:
  sp_group <- c(
    ints[, .(out = unique(group1_id)), by = sp1_id][order(sp1_id)][['out']],
    ints[, .(out = unique(group2_id)), by = sp2_id][order(sp2_id)][['out']]
  )
  site_type <- ints[, .(out = unique(type_id)), by = site_id][order(site_id)][['out']]
  
  # Compute standardized product of suitabilities:
  SS_mean <- mean(ints[['sp1_suitability']] * ints[['sp2_suitability']])
  SS_sd <- sd(ints[['sp1_suitability']] * ints[['sp2_suitability']])
  ints[, SS := scale(sp1_suitability * sp2_suitability)]
  
  # Return list of data:
  list(
    nb_sites = length(site_names),
    nb_spp = length(sp_names),
    nb_int = nrow(ints),
    nb_types = length(type_names),
    nb_groups = length(group_names),
    Y_array = ints[, .(Y, site_id, sp1_id, sp2_id, N)],
    site_type = site_type,
    sp_group = sp_group,
    SS = ints$SS,
    SS_mean = SS_mean,
    SS_sd = SS_sd,
    site_names = site_names,
    sp_ids = sp_ids,
    sp_names = sp_names,
    type_names = type_names,
    group_names = group_names,
    group_types = group_types
  )
}

prepare_stan_start_values.all_binom_bioclim <- function(
  interactions, min_nb_ints = 100, min_bioclim_occs = NULL, collec = FALSE
) {
  data <- prepare_stan_data.all_binom_bioclim(
    interactions, min_nb_ints, min_bioclim_occs, collec
  )
  list(
    alpha = 0,
    lambda_bar = rep(0, data$nb_types),
    zbeta = rep(0, data$nb_sites),
    zgamma = rep(0, data$nb_spp),
    zlambda = rep(0, data$nb_sites),
    sigma_beta = 0.1,
    sigma_gamma = rep(0.1, 2 * data$nb_types),
    sigma_lambda = rep(0.1, data$nb_types)
  )
}
