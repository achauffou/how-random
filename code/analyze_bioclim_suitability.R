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
    unique(by = "cell")
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
    unique(by = "cell")
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
