# Download files from the internet =============================================
download_from_url <- function(download_url, dest_file, download_date = NA) {
  # Create directory where the file should be stored:
  dir.create(dirname(dest_file), showWarnings = FALSE, recursive = TRUE)
  
  # Download the file:
  download.file(download_url, dest_file, quiet = QUIET_DOWNLOADS)
  
  # Return invisibly the destination file:
  dest_file
}


# Download raw Web of Life data ================================================
#' List networks from a given interaction types available in Web of Life
#' 
#' @description 
#' List all the networks that correspond to the specified interaction type. 
#' The download date specified is not taken into account (changing it ensures 
#' that the data is downloaded again). Possible interaction types are:
#' * 3: plant-ant
#' * 5: pollination
#' * 6: seed dispersal
#' * 8: host-parasite
#' * 7: food webs
#' * 10: plant-herbivore
#' * 11: anemone-fish
#' * All
#' 
download_wol_networks_list <- function(
  interaction_type = "All", download_date = NA
) {
  # Resolve interaction type name:
  interaction_id <- interaction_type
  if (is.character(interaction_type)) {
    if (interaction_type != "All") {
      interaction_id <- switch(
        interaction_type,
        "plant-ant" = 3,
        "pollination" = 5,
        "seed-dispersal" = 6,
        "host-parasite" = 8,
        "food-webs" = 7,
        "plant-herbivore" = 10,
        "anemone-fish" = 11,
        "All" = "All"
      )
    }
  }
  
  # Get list of data and parse it with JSON:
  networks_list <- "http://www.web-of-life.es/networkslist.php?type=" %>%
    paste0(interaction_id, "&data=All") %>%
    readLines(warn = FALSE) %>%
    glue::glue_collapse() %>%
    rjson::fromJSON()
  
  # Clean network list to remove networks with no species (root networks):
  networks_list[sapply(networks_list, function (x) x$root) == 0]
}

#' Download Web of Life networks raw data
#' 
download_wol_networks_raw_archive <- function(
  networks_list, dest_file, format = "csv", with_species_names = TRUE
) {
  # Prepare download url:
  networks <- lapply(networks_list, function(x) x[['networkName']]) %>% 
    unlist() %>% 
    paste(collapse = ",")
  species <- FALSE
  if (with_species_names) {species <- "yes"}
  
  # Download data:
  "http://www.web-of-life.es/map_download_fast2.php?" %>%
    paste0("format=", format) %>%
    paste0("&networks=", networks) %>%
    paste0("&species=", species) %>%
    paste0("&type=&data=&speciesrange=&interactionsrange=") %>% 
    paste0("&searchbox=&checked=") %>%
    download_from_url(dest_file, download_date)
}


# Download GBIF occurrence data ================================================
#' Get a data.table with the names to propose and download on GBIF
#' 
select_names_to_gbif_suggest <- function(species, taxonomic_dict, 
                                         min_locs, aggregation_level) {
  # Verified species to download (with nb_locations > min_locations):
  verified_names <- species[
    , .(final_name, loc_id, kingdom), by = .(final_name, loc_id, kingdom)
  ][
    , .(nb_locs = .N), by = .(final_name, kingdom)
  ][
    nb_locs >= min_locs,
  ][, .(final_name, kingdom)]
  
  # Names to suggest:
  names_to_suggest <- taxonomic_dict[, ':='(
    final_name = switch(aggregation_level, "species" = verified_species, 
                         "genus" = verified_genus, verified_name),
    kingdom = verified_kingdom
  )] %>% merge(verified_names, by = c("final_name", "kingdom"))
  names_to_suggest[, .(
    verified_name = final_name, 
    kingdom,
    proposed_name
  )] %>% unique()
}
