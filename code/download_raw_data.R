# Download files from the internet =============================================
download_from_url <- function(download_url, dest_file, download_date = NA) {
  # Create directory where the file should be stored:
  dir.create(dirname(dest_file), showWarnings = FALSE, recursive = TRUE)
  
  # Skip download if the file has been downloaded/created recently:
  if (file.exists(dest_file) & !is.na(download_date)) {
    if (as.Date(file.info(dest_file)$ctime) >= download_date) {
      return(dest_file)
    }
  }
  
  # Download the file if necessary:
  download.file(download_url, dest_file)
  
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
  networks_list, dest_file, format = "csv", with_species_names = TRUE,
  download_date = NA
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


# Download rnaturalearth land data =============================================
#' Download rnaturalearth land data and save it as an R object
#' 
download_rnaturalearth_land_data <- function(
  dest_folder, scale = 10, download_date = NA
) {
  # Output file name:
  dest_file <- dest_folder %>%
    file.path(paste0("rnaturalearth_land_data_", scale, ".rds"))
  
  # Skip download if the file has been downloaded/created recently:
  if (file.exists(dest_file) & !is.na(download_date)) {
    if (as.Date(file.info(dest_file)$ctime) >= download_date) {
      return(dest_file)
    }
  }
  
  # Download data:
  data <- rnaturalearth::ne_download(
    type = "land", category = "physical", returnclass = "sp", scale = scale
  )
  
  # Save data as an R object:
  saveRDS(data, file = dest_file)
  dest_file
}


# Get GBIF keys to downloads ===================================================
#' Get a data.table with the names to propose and download on GBIF
#' 
select_names_to_gbif_suggest <- function(species, taxonomic_dict, 
                                         min_locs, aggregation_level) {
  # Verified species to download (with nb_locations > min_locations):
  verified_names <- species[
    , .(final_id, final_name, loc_id, kingdom), 
    by = .(final_id, final_name, loc_id, kingdom)
  ][
    , .(nb_locs = .N), by = .(final_id, final_name, kingdom)
  ][
    nb_locs >= min_locs,
  ][, .(final_id, final_name, kingdom)]
  
  # Names to suggest:
  names_to_suggest <- taxonomic_dict[, ':='(
    final_name = switch(aggregation_level, "species" = verified_species, 
                         "genus" = verified_genus, verified_name),
    kingdom = verified_kingdom
  )] %>% merge(verified_names, by = c("final_name", "kingdom"))
  names_to_suggest[, .(
    verified_id = final_id,
    verified_name = final_name,
    proposed_name = remove_abbreviations(proposed_name),
    proposed_rank = verified_rank,
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
    nb_cores <- get_nb_cpus()
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
select_gbif_keys_to_download <- function(
  gbif_names_to_suggest, 
  gbif_keys_dict,
  accepted_ranks = c("species", "subspecies")
) {
  # Remove unranked keys and keys from very high orders:
  gbif_keys_dict %<>% .[gbif_rank %in% accepted_ranks,]
  
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
      , .(verified_id, verified_name, proposed_name, proposed_rank, 
          proposed_kingdom, gbif_key, gbif_canonical_name, gbif_kingdom, 
          gbif_rank)
    ]
}


# Download raw GBIF occurrences ================================================
#' Download GBIF occurrence data for given keys (check for existing downloads)
#' 
#' Given a vector of GBIF keys to download, this function checks for each key 
#' whether it is already part of a downloaded dataset. The function downloads 
#' GBIF data for the keys that have not already been downloaded and returns the
#' file names of all archives containing occurrences for the given keys.
#' 
get_gbif_occurrences <- function(
  keys, gbif_data_folder, cache_path, min_date = NA
) {
  # Get only unique keys:
  keys <- unique(keys)
  
  # Update downloads cache:
  downloads_cache <- update_gbif_downloads_cache(gbif_data_folder, 
                                                 cache_path, min_date)
  
  # Get new keys to request:
  new_keys <- keys %>%
    setdiff(downloads_cache[!is.na(request_id)][['gbif_key']])
  
  if (length(new_keys) > 0) {
    # Prepare requests for GBIF keys that do not exist:
    new_requests <- new_keys %>%
      prepare_gbif_requests()
    
    # Request keys that do not exist:
    request_gbif_downloads(new_requests, cache_path)
  }
  
  # Download pending occurrences with successful requests:
  download_pending_gbif_occurrences(keys, gbif_data_folder, cache_path)
  
  # Update downloads cache:
  downloads_cache <- update_gbif_downloads_cache(gbif_data_folder, cache_path, min_date)
  
  # Check if some files have not been downloaded:
  keys_not_found <- downloads_cache[gbif_key %in% keys & is.na(download_file)][['gbif_key']]
  if (length(keys_not_found) > 0) {
    warning(paste(
      "Occurrences for GBIF keys", 
      keys_not_found,
      "could not be retrieved!",
      collapse = ", "
    ))
  }
  
  # Return file names containing the occurrences for the desired keys:
  downloads_cache[gbif_key %in% keys & !is.na(download_file)] %>% 
    .[['download_file']] %>% 
    unique() %>%
    file.path(gbif_data_folder, .)
}

#' Empty cache of GBIF downloads
#' 
empty_gbif_downloads_cache <- function() {
  data.table::data.table(
    gbif_key = integer(),
    request_id = character(),
    download_file = character(),
    request_date = character(),
    download_date = character()
  )
}

#' Update GBIF downloads cache based by removing old, outdated and invalid info
#' 
#' This function remove rows from the download cache that do not correspond to
#' the actual download situation. It first removes any download/request that is 
#' older than the manual download date. It also checks that downloaded files 
#' are present in the download directory and that requests that have not been 
#' downloaded can still be.
#' 
update_gbif_downloads_cache <- function(gbif_data_folder, cache_path, min_date = NA) {
  # Open cache file for GBIF downloads or create it if it does not exist:
  if(file.exists(cache_path)) {
    cache <- data.table::fread(
      cache_path, na.strings = c("", "NA"), 
      colClasses = list(integer = 1, character = 2:5)
    )
  } else {
    dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
    cache <- empty_gbif_downloads_cache()
  }
  
  # Remove any duplicate:
  cache %<>% unique(by = "gbif_key")
  
  # Remove requests/downloads that are too old:
  cache <- cache[request_date >= min_date,]
  cache[download_date < min_date, ':='(download_file = NA, download_date = NA)]
  
  # Check that download files still exist:
  cache[
    !is.na(download_file),
  ][
    !file.exists(file.path(gbif_data_folder, download_file)), 
    ':='(download_file = NA, download_date = NA)
  ]
  
  # Check that pending downloads can still be downloaded:
  request_IDs_to_download <- cache[
    !is.na(request_id) & is.na(download_file),
  ][['request_id']] %>% unique()
  download_info <- rgbif::occ_download_list() %$% 
    results %>% 
    data.table::as.data.table()
  downloadable_request_IDs <- download_info[
    key %in% request_IDs_to_download,][as.Date(eraseAfter) > Sys.Date()]
  outdated_keys <- setdiff(request_IDs_to_download, downloadable_request_IDs)
  cache <- cache[!request_id %in% outdated_keys]
  
  # Save updated cache:
  data.table::fwrite(cache, cache_path)
  cache
}

#' Prepare requests to download GBIF data and return keys and requests
#' 
prepare_gbif_requests <- function(keys) {
  # Function to construct GBIF request:
  construct_gbif_request <- function(keys) {
    rgbif::occ_download_prep(
      rgbif::pred_in("taxonKey", keys),
      rgbif::pred("hasCoordinate", TRUE),
      rgbif::pred("hasGeospatialIssue", FALSE),
      rgbif::pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "OBSERVATION", 
                                        "PRESERVED_SPECIMEN")),
      rgbif::pred_gt("year", 1945)
    )
  }
  
  # Construct initial request:
  request <- construct_gbif_request(keys)
  request_length <- request %$%
    request %>%
    rgbif:::check_inputs() %>%
    as.character() %>%
    nchar()
  
  # If initial request is too long, split it:
  if (request_length > 12000) {
    n_chunks <- request_length %>%
      divide_by(12000) %>%
      ceiling() %>%
      add(1)
    split_keys <- split(keys, cut(seq_along(keys), n_chunks)) %>% unname()
    requests <- lapply(split_keys, construct_gbif_request)
    return(list(split_keys, requests))
  } else {
    return(list(list(keys), list(request)))
  }
}

#' Request prepared GBIF downloads
#' 
request_gbif_downloads <- function(requests, cache_path) {
  # Request prepared GBIF downloads (cancel requests on exit):
  tryCatch({
    message(paste("Requesting", length(requests[[2]]), "downloads from GBIF..."))
    request_ids <- rgbif::occ_download_queue(.list = requests[[2]])
  }, finally = {suppressMessages({rgbif::occ_download_cancel_staged()})})
  
  # If requests were successful, add their IDs and date to downloads cache:
  if (exists("request_ids")) {
    if (length(request_ids) == length(requests[[1]])) {
      # Function to make new cache rows for a single download:
      make_new_cache_rows <- function(nb, requests, request_ids) {
        data.table::data.table(
          gbif_key = requests[[1]][[nb]],
          request_id = as.character(request_ids[[nb]]),
          download_file = NA_character_,
          request_date = Sys.Date(),
          download_date = NA_character_
        )
      }
      
      # Prepare new downloads cache rows:
      1:length(requests[[1]]) %>%
        lapply(make_new_cache_rows, requests, request_ids) %>% 
        data.table::rbindlist() %>% 
        data.table::fwrite(cache_path, append = TRUE)
      return()
    }
  }
  warning("Not all GBIF requests were successful!")
}

#' Download pending requests from GBIF
#' 
download_pending_gbif_occurrences <- function(keys, gbif_data_folder, cache_path) {
  # Read cache:
  cache <- data.table::fread(
    cache_path, na.strings = c("", "NA"), 
    colClasses = list(integer = 1, character = 2:5)
  )
  
  # Get request IDs to download:
  new_downloads <- cache[gbif_key %in% keys & is.na(download_file), 
                         .(request_id)] %>% unique()
  
  # Find download file for each download:
  download_info <- rgbif::occ_download_list() %$% 
    results %>% 
    data.table::as.data.table()
  new_downloads <- download_info[, .(
    request_id = key, 
    new_download_link = downloadLink, 
    new_download_file = basename(downloadLink),
    new_download_date = Sys.Date()
  )] %>% 
    merge(new_downloads, by = "request_id", all.y = TRUE)
  
  # Download pending requests:
  download_pending_request <- function(url, file) {
    success <- FALSE
    trial <- 0
    message(paste(
      "Downloading GBIF request to", file.path(gbif_data_folder, file), "..."
    ))
    while(success == FALSE && trial < 10) {
      tryCatch({
        download_from_url(url, file.path(gbif_data_folder, file))
      }, error = function(e) {}, warning = function(w) {}, finally = {})
      if (file.exists(file.path(gbif_data_folder, file))) success <- TRUE
      Sys.sleep(0.1)
      trial <- trial + 1
    }
    if (success == TRUE) {
      return(file)
    } else {
      return(NA)
    }
  }
  successful_downloads <- purrr::map2_chr(
    new_downloads[!is.na(new_download_link)][['new_download_link']],
    new_downloads[!is.na(new_download_link)][['new_download_file']],
    ~ download_pending_request(.x, .y)
  )
  
  # Update downloads cache:
  new_cache <- new_downloads[new_download_file %in% successful_downloads] %>%
    merge(cache, by = "request_id", all.y = TRUE)
  new_cache[, ':='(
    download_file = ifelse(is.na(download_file), new_download_file, 
                           download_file),
    download_date = ifelse(is.na(download_date), 
                           as.character(new_download_date), download_date)
  )]
  new_cache[
    , .(gbif_key, request_id, download_file, request_date, download_date)
  ] %>% data.table::fwrite(cache_path)
}
