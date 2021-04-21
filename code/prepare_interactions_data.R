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
    is_identified = TRUE,
    is_too_long = FALSE,
    needs_distinction = FALSE,
    distinction_string = NA_character_
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
}

#' Remove abbraviations from species name
#' 
remove_abbreviations <- function(sp_name) {
  sp_name %>%
    stringr::str_replace(stringr::regex("\\s*var\\."), "") %>%
    stringr::str_replace(stringr::regex("\\s*aff\\."), "") %>%
    stringr::str_replace(stringr::regex("\\s*cf\\."), "") %>%
    stringr::str_replace(stringr::regex("\\s*subsp\\."), "") %>%
    stringr::str_replace(" x .+", "")
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
