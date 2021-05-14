# Prepare interactions metadata ================================================
#' Correct manually some locations and add location ID to Web of Life metadata
#' 
create_wol_metadata_loc_id <- function(metadata, manual_loc) {
  # Manual locations corrections:
  metadata <- manual_loc[, .(net_name, man_lat = lat, man_lon = lon)] %>%
    merge(metadata, by = "net_name", all.y = TRUE)
  metadata[!is.na(man_lat), lat := man_lat]
  metadata[!is.na(man_lon), lat := man_lon]
  metadata[, ':='(man_lat = NULL, man_lon = NULL)]
  
  # Add column with location ID:
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
    stringr::str_detect(proposed_name, "^(Unidentified|Undefined|Unientified)"),
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
    stringr::str_replace(stringr::regex("\\s*ssp\\."), "") %>%
    stringr::str_replace(stringr::regex("\\s*spp\\."), "") %>%
    stringr::str_replace(stringr::regex("\\s*sp\\."), "") %>%
    stringr::str_replace(stringr::regex("\\s*indet\\."), "") %>%
    stringr::str_replace(stringr::regex("^ "), "") %>%
    stringr::str_replace(stringr::regex(" $"), "") %>%
    stringr::str_replace(" x .+", "") %>%
    stringr::str_replace(" X ", " ")
}

#' Merge proposed taxonomic names from Web of Life with verified information
#' 
get_verified_names <- function(proposed_species, dict) {
  # Join proposed names with taxonomic dictionary:
  proposed_species %>%
    merge(dict, all.x = TRUE, by = "proposed_name")
}

#' Select and aggregate the most plausible taxonomic data for all species
#' 
#' For all species entries, select the most plausible entry from the verified 
#' names that matches the raw species name. To decide which entry is most 
#' plausible, the function follows the given hierarchy of plausible kingdoms 
#' for the given functional type. Plausible kingdoms must be passed as a named 
#' list. Additionally, entry name is aggregated to the desired taxonomic level 
#' (genus, species or subspecies). Finally, additional columns for final name 
#' and selection flags are added.
#' 
select_verified_species <- function(
  species, 
  verified_names, 
  plausible_kingdoms,
  aggregation_level = "species"
) {
  # Function to get the priority of a verified kingom:
  get_priority <- function(kingdom, fun_group) {
    if (is.na(kingdom)) return(20)
    res <- match(kingdom, plausible_kingdoms[[fun_group]])
    if (is.na(res)) return(10)
    return(res)
  }
  
  # Join the species entry with verified entries matching the raw name:
  species %<>% merge(verified_names, by = "raw_name", all.x = TRUE)
  species %<>% .[, .(
    int_type,
    net_name,
    raw_name,
    final_id = switch(aggregation_level, "species" = sp_id, 
                      "genus" = genus_id, taxon_id),
    final_name = switch(aggregation_level, "species" = verified_species, 
                        "genus" = verified_genus, verified_name),
    kingdom = verified_kingdom,
    loc_id,
    fun_group,
    verified_name,
    is_verified, 
    verification_status
  )]
  
  # Set priority of the entry:
  species[, priority := purrr::map2(kingdom, fun_group, get_priority) %>% 
    as.integer()]
  
  # Flag species with duplicate entries:
  species[
    priority < 10,
  ][
    , selection_flag := ifelse(.N > 1, "Several plausible kingdoms", 
                               NA_character_), 
    by = .(net_name, raw_name)
  ]
  
  # Remove species with duplicate entries:
  species %<>% data.table::setkey(priority) %>%
    unique(by = c("net_name", "raw_name", "fun_group"))
  
  # Flag out species with unknown or implausible kingdom:
  species[priority == 10, selection_flag := "Implausible kingdom"]
  species[priority == 20, selection_flag := "Unknown kingdom"]
  species[, priority := NULL]
}

#' Detect networks with species belonging to incompatible functional groups
#' 
detect_problematic_networks <- function(species, metadata) {
  species[
    !is.na(final_name) & !int_type %in% "Food Webs",
  ][
    , .(problem = length(unique(.SD[['fun_group']])) > 1), by = .(net_name, final_name)
  ][
    problem == TRUE,
  ][['net_name']] %>% 
    c(metadata[is.na(lat) | is.na(lon)][['net_name']]) %>%
    unique()
}

#' Remove species with NA as final name or from problematic networks
#' 
remove_problematic_species <- function(species, problematic_networks) {
  all_net_names <- species[['net_name']] %>% unique()
  problematic_networks <- sapply(
    problematic_networks, 
    function(x) stringr::str_extract(all_net_names, paste0("^", x, ".*$"))
  )
  species[
    !net_name %in% problematic_networks,
  ][
    !is.na(final_name),
  ]
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

#' Get interactions data as a data.table
#' 
#' Remove problematic networks, format as data.table, add location ID,
#' assign cleaned species names, remove species with implausible kingdom and
#' undefined species, and merge duplicate observations at the same location.
#' 
get_wol_interactions <- function(networks, metadata, species, fun_groups, 
                                 problematic_networks = NA) {
  # Remove problematic networks:
  problematic_networks <- sapply(
    problematic_networks, 
    function(x) stringr::str_extract(names(networks), paste0("^", x, ".*$"))
  )
  networks %<>% .[!names(.) %in% problematic_networks]
  
  # Function to format matrix as data.table:
  interactions_as_dt <- function(x) {
    x %>%
      as.table() %>%
      data.table::as.data.table() %>%
      data.table::setnames(c("sp1_raw_name", "sp2_raw_name", "int_weight")) %>%
      .[, ':='(
        sp1_raw_name = as.character(sp1_raw_name),
        sp2_raw_name = as.character(sp2_raw_name),
        int_weight = as.numeric(as.character(int_weight))
      )]
  }
  
  # Remove implausible kingdoms from species:
  species %<>% .[!selection_flag %in% "Implausible kingdom", 
                 .(net_name, raw_name, final_id, final_name, fun_group, kingdom)]
  
  # Prepare interactions data.table:
  interactions <- networks %>%
    # Remove empty networks:
    purrr::discard(~length(.) == 0) %>%
    # Format networks as data.table:
    purrr::map_df(interactions_as_dt, .id = "net_name") %>%
    # Add metadata about location ID and interaction type:
    merge(metadata[, .(net_name, loc_id, int_type)], by = "net_name") %>%
    # Add functional group of each species:
    merge(fun_groups[is_row == TRUE,], by = "int_type") %>%
    .[, ':='(sp1_fun_group = fun_group, fun_group = NULL, is_row = NULL)] %>%
    merge(fun_groups[is_row == FALSE,], by = "int_type") %>%
    .[, ':='(sp2_fun_group = fun_group, fun_group = NULL, is_row = NULL)] %>%
    # Add cleaned species names:
    merge(species[, .(
      net_name, sp1_raw_name = raw_name, sp1_name = final_name, 
      sp1_fun_group = fun_group, sp1_id = final_id
    )], by = c("net_name", "sp1_raw_name", "sp1_fun_group"), all.x = TRUE) %>%
    merge(species[, .(
      net_name, sp2_raw_name = raw_name, sp2_name = final_name, 
      sp2_fun_group = fun_group, sp2_id = final_id
    )], by = c("net_name", "sp2_raw_name", "sp2_fun_group"), all.x = TRUE) %>%
    .[!is.na(sp1_name) & !is.na(sp2_name),]
  
  # Merge duplicate interactions at the same location:
  data.table::setkey(interactions, loc_id, int_type, sp1_fun_group, sp1_id, 
                     sp1_name, sp2_fun_group, sp2_id, sp2_name)
  interactions[
    , .(int_weight = sum(.SD[['int_weight']])), 
    by = .(loc_id, int_type, sp1_fun_group, sp1_id, sp1_name, sp2_fun_group, 
           sp2_id, sp2_name)]
}
