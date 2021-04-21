# Prepare interactions metadata ================================================
#' Add a column with location ID to Web of Life metadata
#' 
create_wol_metadata_loc_id <- function(metadata) {
  metadata[order(net_name), ':='(loc_id = paste0("wol_", .GRP)), by = .(lat, lon)]
}


# Prepare and resolve species names ============================================
#' Return a data.table with raw species names from Web of Life networks
#' 
#' Return a data.table with raw Web of Life species names, and for each record 
#' the name of its network and whether it is a row or column.
#' 
get_raw_wol_species <- function(networks) {
  # Function that returns a data table of species names in network:
  single_network_spp <- function(network) {
    row_spp <- data.table::data.table(
      raw_sp_name = rownames(network), 
      is_row = TRUE
    )
    col_spp <- data.table::data.table(
      raw_sp_name = colnames(network), 
      is_row = FALSE
    )
    rbind(row_spp, col_spp)
  }
  
  # Map over all networks to create the data.table with raw species names:
  networks %<>%
    purrr::discard(~length(.) == 0) %>% # remove any empty network
    purrr::map_dfr(single_network_spp, .id = "net_name")
}

#' Add some metadata to the Web of Life species data.table
#' 
add_metadata_to_wol_species <- function(species, metadata, fun_groups_info) {
  # Add interaction type and location ID to species data.table:
  species %<>% 
    merge(metadata[, .(net_name, int_type, loc_id)], 
          by = c("net_name"), all.x = TRUE)
  
  # Add guild depending on the interaction type:
  species %<>% 
    merge(fun_groups_info, by = c("int_type", "is_row"), all.x = TRUE)
  species[, is_row:=NULL]
}

#' Prepare unchecked species names dictionary
#' 
#' Implement manual species names corrections, remove abbreviations and 
#' aggregate input species names to prepare the cleaning dictionary.
#' 
prepare_unchecked_species_dict <- function(species, manual_corrections = NULL) {
  # Remove duplicate entries:
  species %<>% unique(by = "raw_sp_name")
  
  # Implement manual corrections:
  if (nrow(data.table::as.data.table(manual_corrections)) > 0) {
    dict <- species %>%
      .[, .(raw_sp_name)] %>%
      merge(manual_corrections, by = c("raw_sp_name"), all.x = TRUE)
  } else {
    dict <- species %>%
      .[, .(raw_sp_name)] %>%
      .[, ':='(manual_sp_name = NA_character_, manual_comments = NA_character_)]
  }
  
  # Create all dict columns with default values:
  dict[, ':='(
    prior_sp_name = NA_character_,
    is_valid = TRUE,
    validity_status = NA_character_,
    proposed_sp_name = NA_character_,
    proposed_rank = NA_character_,
    is_identified = TRUE,
    is_too_long = FALSE,
    needs_distinction = FALSE,
    distinction_string = NA_character_,
    proposed_subspecies = NA_character_
  )]
  
  # Create column for proposed species name (the one to query):
  dict[, prior_sp_name := 
         ifelse(is.na(manual_sp_name), raw_sp_name, manual_sp_name)]
  
  # Remove abbreviations from proposed name:
  dict[, proposed_sp_name := remove_abbreviations(prior_sp_name)]
  
  # Detect and flag unidentified species:
  dict[
    stringr::str_detect(proposed_sp_name, "^[Unidentified|Undefined|Unientified]"),
    ':='(
      is_identified = FALSE,
      is_valid = FALSE,
      validity_status = "Unidentified",
      proposed_sp_name = NA_character_
    )
  ]
  
  # Detect and flag entries with several unidentified species of the same genus 
  # in the same network:
  dict[
    stringr::str_detect(proposed_sp_name, "sp[0-9]+.*"),
    ':='(
      distinction_string = proposed_sp_name %>%
        stringr::str_extract("sp[0-9]+.*") %>%
        stringr::str_replace(" ", "_"),
      ambiguous_case = proposed_sp_name %>%
        stringr::str_extract(".*sp[0-9]+.*") %>%
        stringr::str_replace(" ", "_") %>%
        stringr::str_replace("sp[0-9]+", ""),
      proposed_sp_name = proposed_sp_name %>%
        stringr::str_replace(" sp[0-9]+.*", ""),
      needs_distinction = TRUE
    )
  ]
  unique_ambiguous_cases <- dict[, .N, by = .(ambiguous_case)][N == 1, ][['ambiguous_case']]
  dict[
    ambiguous_case %in% unique_ambiguous_cases,
    ':='(
      needs_distinction = FALSE,
      distinction_string = NA_character_
    )
  ][, ambiguous_case := NULL]

  # Aggregate subspecies and indicate a priori rank:
  dict[
    stringr::str_count(proposed_sp_name, " ") > 2,
    ':='(
      is_too_long = TRUE,
      is_valid = FALSE,
      validity_status = "Too long",
      proposed_sp_name = NA_character_
    )
  ][
    stringr::str_count(proposed_sp_name, " ") == 1,
    proposed_rank := 'species'
  ][
    stringr::str_count(proposed_sp_name, " ") == 0,
    proposed_rank := 'genus'
  ][
    stringr::str_count(proposed_sp_name, " ") == 2,
    ':='(
      proposed_rank = 'subspecies',
      proposed_subspecies = proposed_sp_name %>%
        purrr::map(~stringr::str_split(., " ", simplify = TRUE)[3]) %>% 
        unlist(),
      proposed_sp_name = proposed_sp_name %>%
        purrr::map(~stringr::str_split(., " ", simplify = TRUE)[c(1,2)] %>%
                     paste(collapse = " ")) %>% 
        unlist()
    )
  ][is_valid == TRUE, validity_status := "Valid"]
}

#' Remove abbraviations from species name
#' 
remove_abbreviations <- function(sp_name) {
  sp_name %>%
    stringr::str_replace(stringr::regex("\\s*var\\."), "") %>%
    stringr::str_replace(stringr::regex("\\s*aff\\."), "") %>%
    stringr::str_replace(stringr::regex("\\s*cf\\."), "") %>%
    stringr::str_replace(stringr::regex("\\s*subsp\\."), "") %>%
    stringr::str_replace(stringr::regex("\\s*sp\\."), "") %>%
    stringr::str_replace(stringr::regex("\\s*indet\\."), "") %>%
    stringr::str_replace(stringr::regex("^ "), "") %>%
    stringr::str_replace(stringr::regex(" $"), "") %>%
    stringr::str_replace(" x .+", "")
}

#' Automatically check species names on ITIS and GNR (use cache)
#' 
#' Check all proposed species names in a dictionary and return the dictionary 
#' with an additional verified_sp_name column as well as other columns on the 
#' status of verification. To avoid useless verification (involving internet 
#' requests), the function uses the information in the cache file when 
#' the species name has already been assessed.
#' 
check_species_dict <- function(unchecked_dict, synonyms, cache_path) {
  # Read the cache file or create one if it does not exist:
  if(file.exists(cache_path)) {
    cache_dict <- data.table::fread(cache_path)
  } else {
    dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
    cache_dict <- data.table::data.table(
      proposed_sp_name = character(),
      verified_sp_name = character(),
      is_verified = logical(),
      verification_status = character()
    )
  }
  
  # Get proposed names to verify that are not in cache:
  species_names_to_verify <- unchecked_dict[
    !is.na(proposed_sp_name), 
  ][
    !proposed_sp_name %in% cache_dict[['proposed_sp_name']],
  ][['proposed_sp_name']]
  
  # Verify names that are not in cache:
  newly_verified_dict <- species_names_to_verify %>%
    sample(10) %>% # TEMPORARY!!!
    lapply(check_single_species_name, synonyms) %>%
    data.table::rbindlist()
  
  # Append new verifications to cache:
  data.table::fwrite(newly_verified_dict, cache_path, append = TRUE)
  
  # Join verifications to unchecked dictionary:
  full_verification_dict <- rbind(cache_dict, newly_verified_dict)
  unchecked_dict %>%
    merge(full_verification_dict, by = "proposed_sp_name", all.x = TRUE)
}

#' Automatically check a single species name on ITIS and GNR
#' 
check_single_species_name <- function(proposed_sp_name, synonyms) {
  # TEMPORARY!!! Remove once the function is written
  data.table::data.table(
    proposed_sp_name = character(),
    verified_sp_name = character(),
    is_verified = logical(),
    verification_status = character()
  )
}


# Prepare interactions =========================================================
#' Remove supplementary rows/columns from Web of Life networks
#' 
#' Remove any row/column with a name that matches given ones.
#' 
remove_supp_data_from_wol_networks <- function(networks, supp_names) {
  # Remove columns and rows with names that match supplementary data:
  networks %<>%
    purrr::map(~.[!rownames(.) %in% supp_names, !colnames(.) %in% supp_names])
}
