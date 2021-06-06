# Get bioclimatic occurrences at GBIF and Web of Life locations ================
#' Get the GBIF keys corresponding to a species
#' 
get_sp_gbif_keys <- function(name, kingdom, gbif_keys) {
  gbif_keys[verified_name %in% name & proposed_kingdom %in% kingdom][['gbif_key']] %>%
    unique() %>%
    sort()
}

#' Get the bioclimatic conditions at GBIF occurrences of given keys
#' 
get_gbif_bioclim <- function(
  keys, db_path, table_name = "thinned", aggregation_level = "species"
) {
  # Open connection to the database:
  db <- RSQLite::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(RSQLite::dbDisconnect(db))
  
  # Get bioclimatic conditions at occurrences:
  if (aggregation_level %in% "genus") {
    query <- paste0(keys, collapse = ", ") %>% paste0(
      "SELECT DISTINCT * FROM ", table_name, " WHERE ",
      "genusKey IN (", ., ") OR ",
      "speciesKey IN (", ., ") OR ",
      "taxonKey IN (", ., ")", collapse = ", "
    )
  } else if (aggregation_level %in% "species") {
    query <- paste0(keys, collapse = ", ") %>% paste0(
      "SELECT DISTINCT * FROM ", table_name, " WHERE ",
      "speciesKey IN (", ., ") OR ",
      "taxonKey IN (", ., ")", collapse = ", "
    )
  } else {
    query <- paste0(keys, collapse = ", ") %>% paste0(
      "SELECT DISTINCT * FROM ", table_name, " WHERE ",
      "taxonKey IN (", ., ")", collapse = ", "
    )
  }
  RSQLite::dbGetQuery(db, query) %>% 
    data.table::as.data.table() %>% 
    .[, -c("genusKey", "speciesKey", "taxonKey")] %>%
    unique(by = "cell") %>%
    .[complete.cases(.),]
}

#' Get the GBIF bioclimatic conditions of a given species
#' 
get_sp_gbif_bioclim <- function(
  name, kingdom, gbif_keys, db_path, table_name = "thinned", 
  aggregation_level = "species"
) {
  get_sp_gbif_keys(name, kingdom, gbif_keys) %>%
    get_gbif_bioclim(db_path, table_name = "thinned", aggregation_level = "species")
}

#' Get the Web of Life locations of a given species
#' 
get_sp_wol_loc_ids <- function(name, kingdom, wol_species) {
  wol_species[final_name %in% name & kingdom %in% kingdom][['loc_id']] %>%
    unique() %>%
    sort()
}

#' Get the bioclimatic conditions at Web of Life occurrences of a given species
#' 
get_sp_wol_bioclim <- function(name, kingdom, wol_species, wol_bioclim) {
  # Get locations of species:
  wol_bioclim[loc_id %in% get_sp_wol_loc_ids(name, kingdom, wol_species)] %>% 
    .[, -c("loc_id")] %>%
    unique(by = "cell") %>%
    .[complete.cases(.),]
}

#' Get the bioclimatic conditions at all occurrences of a given species
#' 
get_sp_bioclim <- function(
  name, kingdom, wol_species, wol_bioclim, gbif_keys, db_path, 
  table_name = "thinned", aggregation_level = "species"
) {
  rbind(
    get_sp_gbif_bioclim(name, kingdom, gbif_keys, db_path, table_name, aggregation_level),
    get_sp_wol_bioclim(name, kingdom, wol_species, wol_bioclim)
  ) %>% unique(by = "cell")
}

#' Get bioclimatic conditions of all occurrences of all species
#' 
get_all_spp_bioclim <- function(
  wol_bioclim, gbif_keys, db_path, table_name = "thinned", 
  aggregation_level = "species"
) {
  rbind(
    wol_bioclim[complete.cases(wol_bioclim),][, -c("loc_id")],
    gbif_keys[['gbif_key']] %>%
      unique() %>%
      get_gbif_bioclim(db_path, table_name, aggregation_level)
  ) %>% unique(by = "cell")
}

#' Count the number of bioclimatic occurrences of all species
#' 
count_nb_occs_per_species <- function(
  gbif_keys, wol_species, wol_bioclim, db_path, cache_path, 
  table_name = "thinned", aggregation_level = "species"
) {
  # Open cache file or create one if it does not exist:
  if(file.exists(cache_path)) {
    cache <- data.table::fread(
      cache_path, na.strings = c("", "NA"), 
      colClasses = c("character", "character", "integer")
    )
  } else {
    dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
    cache <- data.table::data.table(
      sp_name = character(), 
      sp_kingdom = character(), 
      nb_bioclim_occurrences = integer()
    )
    data.table::fwrite(cache, cache_path)
  }
  
  # Find out which species have not been counted yet:
  spp_to_count <- gbif_keys %>% 
    .[, .(sp_name = verified_name, sp_kingdom = proposed_kingdom)] %>%
    unique() %>%
    data.table::fsetdiff(cache[, .(sp_name, sp_kingdom)])
  
  # Count occurrences of species:
  if (nrow(spp_to_count) > 0) {
    nb_cores <- parallel::detectCores()
    message(paste("Counting all distinct occurrences of", nrow(spp_to_count), "species..."))
    pb <- progress::progress_bar$new(
      total = nrow(spp_to_count),
      format = "  :current/:total :percent |:bar| :elapsedfull",
      incomplete = " ",
      force = TRUE
    )
    parallel::mclapply(1:nrow(spp_to_count), function(x) {
      sp_name <- spp_to_count[x, ][['sp_name']]
      sp_kingdom <- spp_to_count[x, ][['sp_kingdom']]
      nb_occs <- get_sp_bioclim(
        sp_name, sp_kingdom, wol_species, wol_bioclim, gbif_keys, db_path, 
        table_name, aggregation_level
      ) %>% nrow()
      data.table::data.table(sp_name = sp_name, sp_kingdom = sp_kingdom, 
                             nb_bioclim_occurrences = nb_occs) %>%
        data.table::fwrite(cache_path, append = TRUE)
      if(x %% nb_cores == 0) pb$update(x/nrow(spp_to_count))
    }, mc.cores = nb_cores)
    pb$terminate()
    message()
  }
  
  # Return table with occurrences counts for all species:
  data.table::fread(
    cache_path, na.strings = c("", "NA"), 
    colClasses = c("character", "character", "integer")
  )
}


# Calculate niche space and bioclimatic suitability from occurrences ===========
#' Calculate bioclimatic suitability of given occurrences
#' 
#' This function calculates the bioclimatic suitability of a species at given 
#' locations. It takes as input the species and target bioclimatic conditions, 
#' the niche space to use (calculates it based on species bioclimatic 
#' conditions if it is null), as well as the niche space resolution to use. 
#' Species and target bioclimatic occurrences must contain only columns for the 
#' bioclimatic variables and an additional column for the cell number. The 
#' niche size (size of niche containing 90% of occurrences) is also returned. 
#' 
calc_suitability <- function(
  sp_bioclim, target_bioclim, niche_space = NULL, grid_resolution = 200
) {
  # If no niche space was given, calculate it based on the species occurrences:
  if(is.null(niche_space)){
    niche_space <- calc_niche_space(sp_bioclim)
  }
  
  # Get and smooth the species niche:
  sp_niche <- project_bioclim_to_niche_space(sp_bioclim, niche_space) %>%
    realised_niche_density(grid_resolution)
  
  # Calculate niche size:
  area_cell <- raster::res(sp_niche$z) %>% prod()
  niche_size <- sum(as.vector(sp_niche$z > 0.1)) * area_cell
  
  # Get targets bioclimatic conditions as coordinates in the niche space:
  target_niche_coords <- target_bioclim %>%
    project_bioclim_to_niche_space(niche_space) %$%
    lisup
  
  # Calculate the suitability
  suitability <- raster::cellFromXY(sp_niche$z, target_niche_coords) %>%
    raster::extract(sp_niche$z, .)
  
  # Return a data.table with cell, suitability and 90% niche size:
  data.table::data.table(
    cell = target_bioclim[['cell']], 
    suitability = suitability,
    niche_size = niche_size
  )
}

#' Performs a PCA to calculate a bioclimatic niche space from occurrences
#' 
calc_niche_space <- function(bioclim) {
  bioclim %>%
    .[complete.cases(.),] %>%
    .[,-c("cell")] %>%
    ade4::dudi.pca(center = T, scale = T, scannf = F, nf = 2)
}

#' Project bioclimatic conditions into a bidimensional niche space
#' 
project_bioclim_to_niche_space <- function(bioclim, niche_space){
  bioclim %>%
    .[complete.cases(.), ] %>%
    .[,-c("cell")] %>%
    ade4::suprow(niche_space, .)
}

#' Create a grid of occurrences density in a realised niche
#' 
realised_niche_density <- function(realised_niche, grid_resolution){
  realised_niche %$%
    ecospat::ecospat.grid.clim.dyn(lisup, lisup, lisup, R = grid_resolution) %$%
    list(z = z.uncor, w = w)
}


# Analyze the sensitivity of bioclimatic suitability ===========================
#' Compute the climatic suitability of a species for various sample size
#' 
sample_bioclim_suitability_sensitivity <- function(
  sp_name, sp_kingdom, wol_species, wol_bioclim, gbif_keys, db_path, 
  table_name = "thinned", aggregation_level = "species", nb_samples = 10,
  grid_resolution = 200
) {
  # Get species GBIF and Web of Life occurrences:
  gbif_bioclim <- get_sp_gbif_bioclim(
    sp_name, sp_kingdom, gbif_keys, db_path, 
    table_name, aggregation_level
  )
  target_bioclim <- get_sp_wol_bioclim(
    sp_name, sp_kingdom, wol_species, wol_bioclim
  )
  
  # Get sample sizes:
  sampling_levels <- logspace(5, nrow(gbif_bioclim), nb_samples) %>% round()
  
  # Perform sensitivity analysis using the species individual niche:
  nb_cores <- parallel::detectCores()
  message("Performing sensitivity analysis based on individual species niche space...")
  indiv_analysis <- parallel::mclapply(sampling_levels, function(sample_size) {
    res <- gbif_bioclim[sample(.N, sample_size)] %>%
      rbind(target_bioclim) %>%
      unique(by = "cell") %>%
      calc_suitability(target_bioclim, niche_space = NULL,
                       grid_resolution = grid_resolution)
    res[, sample_size := sample_size]
    res
  }, mc.cores = nb_cores) %>% data.table::rbindlist()
  
  # Perform sensitivity analysis using the collective niche of all species:
  message("Performing sensitivity analysis based on collective niche space...")
  collective_niche <- get_all_spp_bioclim(
    wol_bioclim, gbif_keys, db_path, table_name, aggregation_level
  ) %>% calc_niche_space()
  collec_analysis <- parallel::mclapply(sampling_levels, function(sample_size) {
    res <- gbif_bioclim[sample(.N, sample_size)] %>%
      rbind(target_bioclim) %>%
      unique(by = "cell") %>%
      calc_suitability(target_bioclim, niche_space = collective_niche, 
                       grid_resolution = grid_resolution)
    res[, sample_size := sample_size]
    res
  }, mc.cores = nb_cores) %>% data.table::rbindlist()
  
  # Return outcome of sensitivity analyses:
  list(indiv = indiv_analysis, collec = collec_analysis)
}

#' Return a log-spaced sequence
#' 
logspace <- function(x, y, n) exp(seq(log(x), log(y), length.out = n))

#' Calculate the mean absolute error of bioclimatic suitability samples
#' 
calc_bioclim_suitability_sensitivity_errors <- function(sensitivity_samples) {
  # Get the mean absolute error:
  sensitivity_samples[
    , baseline_suitability := .SD[order(-sample_size)][['suitability']][1], 
    by = .(cell)
  ]
  sensitivity_samples[
    , abs_error := abs(.SD[['suitability']] - baseline_suitability), 
    by = .(cell, sample_size)
  ]
  sensitivity_errors <- sensitivity_samples[
    , .(mae = mean(abs_error)), by = .(sample_size)
  ]
  sensitivity_errors[
    , log_prop_occ := log(sample_size / max(.SD[['sample_size']]))
  ]
}

#' Determine the number of occurrences to meet suitability precision criterion
#' 
calc_bioclim_suitability_min_occurrences <- function(
  sensitivity_errors, max_suitability_error
) {
  # Make a simple regression between log_prop_occ and mae:
  threshold <- sensitivity_errors %>% 
    mgcv::gam(mae ~ s(log_prop_occ), data = . , family = "gaussian") %>% 
    broom::augment(
      type.predict = "response", newdata = tibble::tibble(
        log_prop_occ = seq(
          min(sensitivity_errors[['log_prop_occ']]), 
          max(sensitivity_errors[['log_prop_occ']]), 
          length.out = 1000)
      )
    ) %$%
    approx(x = .fitted, y = log_prop_occ, max_suitability_error) %$%
    exp(y) * max(sensitivity_errors[['sample_size']])
  round(threshold)
}
