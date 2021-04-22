# Prepare interactions metadata ================================================
#' Add a column with location ID to Web of Life metadata
#' 
create_wol_metadata_loc_id <- function(metadata) {
  metadata[order(net_name), ':='(loc_id = paste0("wol_", .GRP)), by = .(lat, lon)]
}


# Prepare Web of Life species ==================================================
#' Return a data.table with raw species names from Web of Life networks
#' 
#' Return a data.table with raw Web of Life species names, and for each record 
#' the name of its network and whether it is a row or column.
#' 
get_raw_wol_species <- function(networks) {
  # Function that returns a data table of species names in network:
  single_network_spp <- function(network) {
    row_spp <- data.table::data.table(
      raw_name = rownames(network), 
      is_row = TRUE
    )
    col_spp <- data.table::data.table(
      raw_name = colnames(network), 
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

#' Prepare taxonomic names to propose for taxonomic verification
#' 
#' Implement manual names corrections, remove abbreviations, flag invalid 
#' names and remove distinction string of unidentified species from identified 
#' genera. The function returns a data.table with a proposed_name column that 
#' contains the names to use for taxonomic verification (set to NA for invalid 
#' names).
#' 
prepare_names_to_verify <- function(species, manual_corrections = NULL) {
  # Remove duplicate entries:
  species %<>% unique(by = "raw_name")
  
  # Implement manual corrections:
  if (nrow(data.table::as.data.table(manual_corrections)) > 0) {
    dict <- species %>%
      .[, .(raw_name)] %>%
      merge(manual_corrections, by = c("raw_name"), all.x = TRUE)
  } else {
    dict <- species %>%
      .[, .(raw_name)] %>%
      .[, ':='(manual_name = NA_character_, manual_comments = NA_character_)]
  }
  
  # Create all dict columns with default values:
  dict[, ':='(
    prior_name = NA_character_,
    proposed_name = NA_character_,
    proposed_level = NA_character_,
    is_valid = TRUE,
    validity_status = NA_character_,
    is_identified = TRUE,
    is_too_long = FALSE,
    is_corrected = FALSE,
    had_abbreviations = FALSE
  )]
  
  # Create column for proposed species name (the one to query):
  dict[
    , prior_name := ifelse(is.na(manual_name), raw_name, manual_name)
  ][
    , is_corrected := raw_name != prior_name
  ]
  
  # Remove abbreviations from proposed name:
  dict[
    , proposed_name := remove_abbreviations(prior_name)
  ][
    , had_abbreviations := prior_name != proposed_name
  ]
  
  # Detect and flag unidentified species:
  dict[
    stringr::str_detect(proposed_name, "^[Unidentified|Undefined|Unientified]"),
    ':='(
      is_identified = FALSE,
      is_valid = FALSE,
      validity_status = "Unidentified",
      proposed_name = NA_character_
    )
  ]
  
  # Remove distinction string for unidentified species of identified genus:
  dict[
    , proposed_name := stringr::str_replace(proposed_name, " sp[0-9]+.*", "")
  ]
  
  # Detect and flag names that are too long:
  dict[
    stringr::str_count(proposed_name, " ") > 2,
    ':='(
      is_too_long = TRUE,
      is_valid = FALSE,
      validity_status = "Too long",
      proposed_name = NA_character_
    )
  ]
  
  # Find the hierarchy level of the valid proposed name:
  dict[
    , proposed_level := c('higher', 'species', 'subspecies') %>%
      .[stringr::str_count(proposed_name, " ") + 1]
  ][
    is_valid == TRUE, validity_status := "Valid"
  ]
}

#' Remove abbraviations from taxon name
#' 
remove_abbreviations <- function(name) {
  name %>%
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
