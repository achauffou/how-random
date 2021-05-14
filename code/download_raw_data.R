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


# Get GBIF keys to downloads ===================================================
#' Get a data.table with the names to propose and download on GBIF
#' 
select_names_to_gbif_suggest <- function(species, taxonomic_dict, 
                                         min_locs, aggregation_level) {
  # Verified species to download (with nb_locations > min_locations):
  verified_names <- species[
    , .(final_name, final_rank, loc_id, kingdom), 
    by = .(final_name, loc_id, kingdom)
  ][
    , .(final_rank, nb_locs = .N), by = .(final_name, kingdom)
  ][
    nb_locs >= min_locs,
  ][, .(final_name, final_rank, kingdom)]
  
  # Names to suggest:
  names_to_suggest <- taxonomic_dict[, ':='(
    final_name = switch(aggregation_level, "species" = verified_species, 
                         "genus" = verified_genus, verified_name),
    kingdom = verified_kingdom
  )] %>% merge(verified_names, by = c("final_name", "kingdom"))
  names_to_suggest[, .(
    verified_name = final_name,
    proposed_name = remove_abbreviations(proposed_name),
    proposed_rank = final_rank,
    proposed_kingdom = kingdom
  )] %>% unique()
}

#' Check suggested names in GBIF and returns the corresponding keys (with cache)
#' 
suggest_gbif_names <- function(
  suggested_names, 
  cache_path
) {
  # Create cache file for GBIF keys if it does not exist:
  if(file.exists(cache_path)) {
    cache_keys <- data.table::fread(
      cache_path, na.strings = c("", "NA"), colClasses = list(
        character = c(1, 2, 3, 5, 6, 7),
        integer = 4
      )
    )
  } else {
    dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
    cache_keys <- empty_gbif_keys()
    data.table::fwrite(cache_keys, cache_path)
  }
  
  # Names to suggest:
  names_to_suggest <- suggested_names[
    !is.na(proposed_name),
  ][
    !proposed_name %in% cache_keys[['proposed_name']],
  ] %>% unique()
  
  # Check names that have not been checked:
  if (nrow(names_to_suggest) > 0) {
    nb_cores <- parallel::detectCores()
    message(paste("Checking", nrow(names_to_suggest), "names in GBIF..."))
    pb <- progress::progress_bar$new(
      total = nrow(names_to_suggest),
      format = "  :current/:total :percent |:bar| :elapsedfull",
      incomplete = " ",
      force = TRUE
    )
    new_keys <- 1:nrow(names_to_suggest) %>%
      parallel::mclapply(
        function(x) {
          proposed_name <- names_to_suggest[x][['proposed_name']]
          proposed_rank <- names_to_suggest[x][['proposed_rank']]
          proposed_kingdom <- names_to_suggest[x][['proposed_kingdom']]
          suggest_single_gbif_name(proposed_name, proposed_rank, 
                                   proposed_kingdom, cache_path)
          if(x %% nb_cores == 0) pb$update(x/nrow(names_to_suggest))
        }, mc.cores = nb_cores
      )
    pb$terminate()
    message()
  }
  
  # Return GBIF keys dictionary:
  data.table::fread(
    cache_path, na.strings = c("", "NA"), colClasses = list(
      character = c(1, 2, 3, 5, 6, 7),
      integer = 4
    )
  ) %>% unique() %>% .[!is.na(gbif_key),]
}

#' Create an empty table for GBIF keys
#' 
empty_gbif_keys <- function() {
  data.table::data.table(
    proposed_name = character(),
    proposed_rank = character(),
    proposed_kingdom = character(),
    gbif_key = integer(),
    gbif_canonical_name = character(),
    gbif_kingdom = character(),
    gbif_rank = character()
  )
}

#' Get GBIF keys for a single proposed name
#' 
suggest_single_gbif_name <- function(
  suggested_name, suggested_rank, suggested_kingdom, cache_path
) {
  # Suggest the name to GBIF:
  gbif_results <- rgbif::name_suggest(
    suggested_name, limit = 1000, 
    fields = c("key", "kingdom", "canonicalName", "rank")
  ) %$% data %>% data.table::as.data.table()
  
  # Remove keys from incorrect kingdoms:
  if (nrow(gbif_results) > 0) {
    keys <- switch(suggested_kingdom,
      "Animalia" = gbif_results[kingdom %in% "Animalia",],
      "Plantae" = gbif_results[kingdom %in% "Plantae",],
      "Fungi" = gbif_results[kingdom %in% "Fungi",],
      "Other Eukaryota" = gbif_results[kingdom %in% c("Chromista", "Protozoa"),],
      "Archaea" = gbif_results[kingdom %in% "Archaea",],
      "Bacteria" = gbif_results[kingdom %in% "Bacteria",],
      gbif_results[!kingdom %in% "Viruses",]
    )
    keys[, rank := tolower(rank)]
  } else {keys <- data.table::data.table()}
  
  # Format results:
  if (nrow(keys) > 0) {
    results <- keys[, .(
        proposed_name = suggested_name,
        proposed_rank = suggested_rank,
        proposed_kingdom = suggested_kingdom,
        gbif_key = key,
        gbif_canonical_name = canonicalName,
        gbif_kingdom = kingdom,
        gbif_rank = rank
      )]
  } else {
    results <- data.table::data.table(
      proposed_name = suggested_name,
      proposed_rank = suggested_rank,
      proposed_kingdom = suggested_kingdom,
      gbif_key = NA_integer_,
      gbif_canonical_name = NA_character_,
      gbif_kingdom = NA_character_,
      gbif_rank = NA_character_
    )
  }
  
  # Write results to cache:
  data.table::fwrite(results, cache_path, append = TRUE)
}

#' Select the GBIF keys to download
#' 
select_gbif_keys_to_download <- function(gbif_names_to_suggest, gbif_keys_dict) {
  # Remove unranked keys and keys from very high orders:
  gbif_keys_dict %<>% .[!gbif_rank %in% c("UNRANKED", "OTHER"),]
  
  # Select GBIF keys to include:
  gbif_keys <- gbif_keys_dict[, ':='(
    same_name = (tolower(proposed_name) == tolower(gbif_canonical_name)),
    same_level = stringr::str_count(proposed_name, " ") == 
      stringr::str_count(gbif_canonical_name, " ")
  )][
    , ':='(
      proposed_found = any(.SD[['same_name']] == TRUE),
      only_sub = all(.SD[['same_level']] == FALSE)
    ), by = .(proposed_name)
  ][, ':='(
    necessary_name = (same_name | !proposed_found)
  )][
    (same_level & necessary_name) | only_sub,
  ]
  
  # Merge GBIF keys with verified names:
  gbif_names_to_suggest %>% 
  merge(gbif_keys, by = c("proposed_name", "proposed_rank", "proposed_kingdom")) %>% 
    .[
      , .(verified_name, proposed_name, proposed_rank, proposed_kingdom, gbif_key, 
          gbif_canonical_name, gbif_kingdom, gbif_rank)
    ]
}
  ]
}
