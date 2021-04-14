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
list_wol_networks <- function(interaction_type = "All", download_date = NA) {
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
  "http://www.web-of-life.es/networkslist.php?type=" %>%
    paste0(interaction_id, "&data=All") %>%
    readLines(warn = FALSE) %>%
    glue::glue_collapse() %>%
    rjson::fromJSON()
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
