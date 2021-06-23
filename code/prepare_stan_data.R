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
    fun_args <- list()
  } else {
    fun_args <- spec$data_generation_args
  }
  paste0("prepare_stan_data.", spec$data_generation_model) %>%
    get() %>%
    do.call(args = c(interactions, fun_args)) %>%
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
    fun_args <- list()
  } else {
    fun_args <- spec$data_generation_args
  }
  paste0("prepare_stan_start_values.", spec$data_generation_model) %>%
    get() %>%
    do.call(args = c(interactions, fun_args)) %>%
    saveRDS(file = out_path)
  out_path
}


# Functions to prepare data for specific Stan analyses =========================
prepare_stan_data.pol_logit_03 <- function(
  interactions, min_bioclim_occs = NULL, collec = FALSE
) {
  # Select relevant columns with complete cases:
  if (!is.null(min_bioclim_occs)) {
    ints <- interactions[sp1_nb_bioclim >= min_bioclim_occs & 
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
  
  # Prepare sites, plant and pollinators IDs:
  ints %<>% .[, pla_ID := .GRP, by = .(sp1_id)] %>%
    .[, pol_ID := .GRP, by = .(sp2_id)] %>%
    .[, site_ID := .GRP, by = .(net_id)]
  pla_names <- ints[, .(sp1_name), by = pla_ID][order(pla_ID)][['sp1_name']]
  pla_kingdoms <- ints[, .(sp1_kingdom), by = pla_ID][order(pla_ID)][['sp1_kingdom']]
  pol_names <- ints[, .(sp2_name), by = pol_ID][order(pol_ID)][['sp1_name']]
  pol_kingdoms <- ints[, .(sp2_kingdom), by = pol_ID][order(pol_ID)][['sp1_kingdom']]
  site_names <- ints[, .(net_id), by = site_ID][order(site_ID)][['net_id']]
  
  # Use binary interactions:
  ints[, ':='(Y := ifelse(int_strength > 0, 1, 0))]
  
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
    Y_array = data[, .(Y, site_id, pla_id, pol_id)],
    SS = ints$SS,
    SS_mean = SS_mean,
    SS_sd = SS_sd,
    site_names = site_names,
    pla_names = pla_names,
    pla_kingdoms = pla_kingdoms,
    pol_names = pol_names,
    pol_kingdoms = pol_kingdoms
  )
}

prepare_stan_data.pol_logit_03 <- function(
  interactions, min_bioclim_occs = NULL, collec = FALSE
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
