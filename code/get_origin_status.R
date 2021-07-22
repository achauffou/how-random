# Find out the country and region of each site =================================
#' Get the country or region of each Web of Life site
#' 
get_sites_regions_codes <- function(metadata, wab_countries, wgsrpd_l3, manual_codes) {
  # Determine country/region code using spatial polygons:
  site_codes <- metadata[
    !is.na(lat)
  ][
    , .(lat = unique(lat), lon = unique(lon)), by = .(loc_id)
  ]
  site_codes[, ':='(
    WAB_code = get_wab_code(lon, lat, wab_countries),
    WGSRPD_code = get_wgsrpd_code(lon, lat, wgsrpd_l3)
  )]
  
  # Implement manual corrections:
  if (nrow(manual_codes) > 0) {
    manual_codes %<>% 
      merge(site_codes[, .(loc_id, lat, lon)], by = "loc_id", all.x = TRUE)
    res <- rbind(
      site_codes[!loc_id %in% manual_codes[['loc_id']]],
      manual_codes[, .(loc_id, lat, lon, WAB_code, WGSRPD_code)]
    )
  } else {
    res <- site_codes
  }
  res
}

#' Get WAB country code of given coordinates
#' 
get_wab_code <- function(lon, lat, wab_countries) {
  data.frame(lon = lon, lat = lat) %>%
    sp::SpatialPoints(proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
    sp::over(wab_countries) %$%
    iso_3166_1_
}

#' Get WGSRPD region code of given coordinates
#' 
get_wgsrpd_code <- function(lon, lat, wgsrpd_l3) {
  data.frame(lon = lon, lat = lat) %>%
    sp::SpatialPoints(proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
    sp::over(wgsrpd_l3) %$%
    LEVEL3_COD
}


# Determine neighbour regions for WAB and WGSRPD polygons =====================
#' Get WAB countries neighbours
#' 
get_wab_neighbours <- function(wab_countries) {
  neighbours <- spdep::poly2nb(wab_countries) %>%
    purrr::map(~ if (all(. == 0)) {NA_character_} else {wab_countries$iso_3166_1_[.]})
  names(neighbours) <- wab_countries$iso_3166_1_
  neighbours
}

#' Get WGSRPD regions neighbours
#' 
get_wgsrpd_neighbours <- function(wgsrpd_l3) {
  neighbours <- spdep::poly2nb(wgsrpd_l3) %>%
    purrr::map(~ if (all(. == 0)) {NA_character_} else {wgsrpd_l3$LEVEL3_COD[.]})
  names(neighbours) <- wgsrpd_l3$LEVEL3_COD
  neighbours
}


# Determine the origin status of species at all locations ======================
#' Get the origin status of all species at all locations
#' 
get_spp_origin_status <- function(
  wol_species, taxonomic_dict, sites_regions_codes, wab_neighbours, 
  wgsrpd_neighbours, cache_path, aggregation_level
) {
  # If there is no cache, create one:
  if (!file.exists(cache_path)) {
    cache <- data.table::data.table(
      sp_name = character(),
      sp_kingdom = character(),
      loc_id = character(),
      is_native = logical(),
      status = character(),
      pow_status = character(),
      iucn_status = character(),
      has_conflict = logical()
    )
    data.table::fwrite(cache, cache_path)
  } else {
    cache <- data.table::fread(
      cache_path, na.strings = c("", "NA"), 
      colClasses = list(character = c(1:3, 5:7), logical = c(4, 8))
    )
  }
  
  # Prepare the data table with all species at all locations:
  origins <- wol_species %>%
    .[, .N, by = .(sp_name = final_name, sp_kingdom = kingdom, loc_id)] %>%
    merge(sites_regions_codes, by = "loc_id") %>%
    .[!(is.na(WAB_code) & is.na(WGSRPD_code)), -c("N")]
  
  # Get species/locations combinations that have not been retrieved yet:
  spplocs_to_check <- data.table::fsetdiff(
    origins[, .(sp_name, sp_kingdom, loc_id)],
    cache[, .(sp_name, sp_kingdom, loc_id)]
  )
  
  # Find out the origin status of each species at all its locations:
  if (nrow(spplocs_to_check) > 0) {
    spp_names <- spplocs_to_check[, .N, by = .(sp_name, sp_kingdom)][['sp_name']]
    spp_kingdoms <- spplocs_to_check[, .N, by = .(sp_name, sp_kingdom)][['sp_kingdom']]
    nb_cores <- get_nb_cpus()
    message(paste("Checking origin status of", length(spp_names), "species..."))
    pb <- progress::progress_bar$new(
      total = length(spp_names),
      format = "  :current/:total :percent |:bar| :elapsedfull",
      incomplete = " ",
      force = TRUE
    )
    parallel::mclapply(1:length(spp_names), function(x) {
      tryCatch({
        locs <- spplocs_to_check[
          sp_name %in% spp_names[x] & sp_kingdom %in% spp_kingdoms[x]
        ][['loc_id']] %>% unique()
        check_single_sp_origin_status(
          spp_names[x], spp_kingdoms[x], locs, taxonomic_dict, 
          sites_regions_codes, wab_neighbours, wgsrpd_neighbours, cache_path,
          aggregation_level, nb_cores)
      }, error = function(e) {
        cat("Species sp_name=", spp_names[x],
            ", sp_kingdom=", spp_kingdoms[x],
            " issued an error!\n\t", conditionMessage(e), "\n")
      }, warning = function(w) cat(
        "Species sp_name=", spp_names[x],
        ", sp_kingdom=", spp_kingdoms[x],
        " issued a warning!\n\t", conditionMessage(w), "\n"), finally = {})
      if(x %% nb_cores == 0) pb$update(x/length(spp_names))
    }, mc.cores = nb_cores)
    pb$terminate()
    message()
    
    # Read the updated cache:
    cache <- data.table::fread(
      cache_path, na.strings = c("", "NA"), 
      colClasses = list(character = c(1:3, 5:7), logical = c(4, 8))
    )
  }
  merge(origins[, .(sp_name, sp_kingdom, loc_id)], cache, 
        by = c("loc_id", "sp_name", "sp_kingdom"), all.x = TRUE)
}

#' Check the origin status of a single species at several locations
#' 
check_single_sp_origin_status <- function(
  sp_name, sp_kingdom, locs, taxonomic_dict, sites_regions_codes, 
  wab_neighbours, wgsrpd_neighbours, cache_path, aggregation_level, nb_cores = 1
) {
  # Get names to propose:
  if (aggregation_level == "genus") {
    proposed_names <- taxonomic_dict[
      verified_genus %in% sp_name & verified_kingdom %in% sp_kingdom
    ][['proposed_name']] %>% unique()
  } else if (aggregation_level == "species") {
    proposed_names <- taxonomic_dict[
      verified_species %in% sp_name & verified_kingdom %in% sp_kingdom
    ][['proposed_name']] %>% unique()
  } else {
    proposed_names <- taxonomic_dict[
      verified_name %in% sp_name & verified_kingdom %in% sp_kingdom
    ][['proposed_name']] %>% unique()
  }
  
  # Get species POW and IUCN distribution data:
  pow_distr <- get_pow_distr_data(proposed_names, taxonomic_dict, nb_cores)
  iucn_distr <- get_iucn_distr_data(proposed_names, nb_cores)
  
  # Apply origin status verification to locations where the species is found:
  sp_origins <- locs %>%
    purrr::map_df(
      apply_origin_status_heuristics, sites_regions_codes, wab_neighbours, 
      wgsrpd_neighbours, pow_distr, iucn_distr
    ) %>%
    data.table::as.data.table()
  sp_origins[, ':='(sp_name = sp_name, sp_kingdom = sp_kingdom)]
  data.table::setcolorder(sp_origins, c(
    "sp_name", "sp_kingdom", "loc_id", "is_native", "status", "pow_status", 
    "iucn_status", "has_conflict"
  ))
  data.table::fwrite(sp_origins, cache_path, append = TRUE)
}


# Retrieve species distribution data ===========================================
#' Get POW distribution data of proposed names
#' 
get_pow_distr_data <- function(proposed_names, taxonomic_dict, nb_cores = 1) {
  fqIds <- taxonomic_dict[
    proposed_name %in% proposed_names & !is.na(pow_fqId)][['pow_fqId']] %>%
    unique()
  if (length(fqIds) > 0) {
    all_distr_data <- lapply(fqIds, get_single_pow_distr, nb_cores)
    distr <- list(
      native = unique(c(unlist(sapply(all_distr_data, function(x) x[["native"]])))),
      introduced = unique(c(unlist(sapply(all_distr_data, function(x) x[["introduced"]]))))
    )
  } else {
    distr <- list(native = character(), introduced = character())
  }
  distr
}

#' Get Plants of the World distribution data of a single fqId
#' 
get_single_pow_distr <- function(fqId, nb_cores = 1) {
  success <- FALSE
  trial <- 0
  res <- NULL
  while (success == FALSE & trial < 10) {
    trial <- trial + 1
    Sys.sleep(0.21 * nb_cores)
    tryCatch({
      res <- taxize::pow_lookup(fqId, include = "distribution")
    }, error = function(e) {}, finally = {})
    if (!is.null(res$meta)) success <- TRUE
  }
  if (is.null(res$meta$distribution)) {
    native <- character()
    introduced <- character()
  } else {
    if (!all(names(res$meta$distribution) %in% 
             c("natives", "introduced", "doubtful", "absent", "extinct"))) {
      stop(paste("POW distribution status unaccounted for among", paste(
        names(res[[1]]$distr), collapse =", "
      )))
    }
    native <- as.character(res$meta$distribution$natives[["tdwgCode"]])
    introduced <- as.character(res$meta$distribution$introduced[["tdwgCode"]])
  }
  list(native = native, introduced = introduced)
}

#' Get IUCN distribution data of proposed names:
#' 
get_iucn_distr_data <- function(proposed_names, nb_cores = 1) {
  all_distr_data <- lapply(unique(proposed_names), get_single_iucn_distr, nb_cores)
  distr <- list(
    native = unique(c(unlist(sapply(all_distr_data, function(x) x[["native"]])))),
    introduced = unique(c(unlist(sapply(all_distr_data, function(x) x[["introduced"]]))))
  )
}

#' Get IUCN Red List distribution data of a single proposed name
#' 
get_single_iucn_distr <- function(name, nb_cores) {
  if (stringr::str_count(name, " ") != 1) return(list(native = NULL, invasive = NULL))
  success <- FALSE
  trial <- 0
  res <- NULL
  while (success == FALSE & trial < 10) {
    trial <- trial + 1
    Sys.sleep(2.01 * nb_cores)
    sink(file = file.path(tempdir(), "iucn_sink.txt"), type = "output")
    tryCatch({suppressMessages({suppressWarnings({
      res <- taxize::iucn_summary(name, distr_detail = TRUE)
    })})}, error = function(e) {}, finally = {})
    sink(type = "output")
    if (!is.na(res)) success <- TRUE
  }
  if (length(names(res[[1]]$distr)) == 0) {
    native <- character()
    introduced <- character()
  } else {
    if (!all(names(res[[1]]$distr) %in% c("Native", "Introduced", "Vagrant", 
                                          "Present - Origin Uncertain",
                                          "Regionally Extinct", 
                                          "Possibly Extinct"))) {
      stop(paste("IUCN distribution status unaccounted for among", paste(
        names(res[[1]]$distr), collapse = ", "
      )))
    }
    native <- c(
      as.character(res[[1]]$distr$Native[["code"]]), 
      as.character(res[[1]]$distr$Vagrant[["code"]])
    )
    introduced <- as.character(res[[1]]$distr$Introduced[["code"]])
  }
  list(native = native, introduced = introduced)
}


# Apply origin status heuristics for a single species at a single locations ====
#' Return the origin status of a given species at some location
#' 
apply_origin_status_heuristics <- function(
  loc, sites_regions_codes, wab_neighbours, wgsrpd_neighbours, pow_distr, 
  iucn_distr
) {
  # Find out country/region codes:
  wgsrpd_code <- sites_regions_codes[loc_id == loc][['WGSRPD_code']]
  wab_code <- sites_regions_codes[loc_id == loc][['WAB_code']]
  
  # Get origin status from both databases:
  pow_status <- deduce_origin_status(wgsrpd_code, pow_distr, wgsrpd_neighbours)
  iucn_status <- deduce_origin_status(wab_code, iucn_distr, wab_neighbours)
  
  # Determine the overall status and identify conflicts:
  has_conflict <- FALSE
  if (pow_status %in% "Native") {
    if (iucn_status %in% "Introduced") {
      is_native <- NA
      status <- "Conflict"
      has_conflict <- TRUE
    } else {
      is_native <- TRUE
      status <- "Native"
    }
  } else if (pow_status %in% "Introduced") {
    if (iucn_status %in% "Native") {
      is_native <- NA
      status <- "Conflict"
      has_conflict <- TRUE
    } else {
      is_native <- FALSE
      status <- "Introduced"
    }
  } else if (iucn_status %in% "Native") {
    is_native <- TRUE
    status <- "Native"
  } else if (iucn_status %in% "Introduced") {
    is_native <- FALSE
    status <- "Introduced"
  } else if (pow_status %in% "Neighbour native") {
    is_native <- TRUE
    status <- "Neighbour native"
  } else if (iucn_status %in% "Neighbour native") {
    is_native <- TRUE
    status <- "Neighbour native"
  } else if (pow_status %in% "Neighbour introduced") {
    is_native <- FALSE
    status <- "Neighbour introduced"
  } else if (iucn_status %in% "Neighbour introduced") {
    is_native <- FALSE
    status <- "Neighbour introduced"
  } else if (pow_status == "Unknown" | iucn_status == "Unknown") {
    is_native <- NA
    status <- "Unknown"
  } else {
    is_native <- NA
    status <- "No data"
  }
  
  # Return formatted data.table:
  data.table::data.table(
    loc_id = loc,
    is_native = is_native,
    status = status,
    pow_status = pow_status,
    iucn_status = iucn_status,
    has_conflict = has_conflict
  )
}

#' Return POW or IUCN origin status
#' 
deduce_origin_status <- function(loc_code, distr, neighbours) {
  # Check whether data is available:
  if (is.na(loc_code)) return("No location")
  if (length(distr$native) == 0 & length(distr$introduced) == 0) return("No data")
  
  # Check if the species is registered in the zone as native/introduced:
  if (loc_code %in% distr$native) return("Native")
  if (loc_code %in% distr$introduced) return("Introduced")
  
  # If the species has no status in the zone, try neighbour zones:
  if (loc_code %in% unlist(neighbours[distr$native])) return("Neighbour native")
  if (loc_code %in% unlist(neighbours[distr$introduced])) return("Neighbour introduced")
  
  # By default, return an unknown status:
  "Unknown"
}
