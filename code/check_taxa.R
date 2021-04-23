# Proceed to taxonomic verification using cached dictionary ====================
#' Check proposed taxonomic names and return taxonomic dictionary
#' 
#' This function verifies the validity of proposed taxonomic names and add 
#' those names as well as their verified counterparts and synonyms to the 
#' taxonomic dictionary. To spare needless internet requests and computations,
#' the function uses a cache file in which is stored all names previously
#' proposed. The function returns the entire updated species dictionary.
#' 
check_proposed_names <- function(
  proposed_names, 
  itis_database, 
  cache_path = file.path(tempdir(), "taxonomic_dictionary_cache.csv")
) {
  # Open cached dictionary or create one if it does not exist:
  if(file.exists(cache_path)) {
    cache_dict <- data.table::fread(cache_path, na.strings = c("", "NA"))
  } else {
    dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
    cache_dict <- empty_taxonomic_dict()
    data.table::fwrite(cache_dict, cache_path)
  }
  
  # Get proposed names that are not in cache:
  names_to_verify <- proposed_names[
    !is.na(proposed_name), 
  ][
    !proposed_name %in% cache_dict[['proposed_name']],
  ][['proposed_name']]
  
  # Verify proposed names that are not already in cache:
  names_to_verify %>%
    sample(10) %>% # TEMPORARY!!!
    lapply(check_single_name, itis_database, cache_path)
  
  # Remove duplicates, set taxon IDs and update last proposition date:
  data.table::fread(cache_path, na.strings = c("", "NA")) %>%
    .[order(-verification_time), ] %>%
    unique(by = c("proposed_name", "verified_kingdom", "verified_name")) %>%
    set_taxon_ID() %>%
    .[proposed_name %in% proposed_names, last_proposed_time := character(Sys.time())]
}

#' Create an empty taxonomic dictionary
#' 
empty_taxonomic_dict <- function() {
  data.table::data.table(
    proposed_name = character(),
    verified_name = character(),
    verified_level = character(),
    verified_rank = character(),
    verified_kingdom  = character(),
    is_verified = logical(),
    verification_status = character(),
    gnr_status = character(),
    gnr_match = character(),
    gnr_score = numeric(),
    gnr_source = character(),
    found_in_itis = character(),
    itis_tsn = integer(),
    itis_accepted_tsn = integer(),
    itis_status = character(),
    itis_reason = character(),
    itis_rank = character(),
    itis_kingdom = character(),
    verification_time = character(),
    last_proposed_time = character()
  )
}

#' Check a single proposed taxonomic name (use and save in cache)
#' 
#' Check the proposed name in the taxonomic directory in cache and verify the 
#' name if it is not already present in the dictionary. If the name was not 
#' already in the taxonomic dictionary, update the cache accordingly after 
#' verifying it.
#' 
check_single_name <- function(name, itis_database, cache_path) {
  # Check whether the proposed name is already in the cache:
  cache_dict <- data.table::fread(cache_path, na.strings = c("", "NA"))
  
  # If the name is not already in cache, verify it:
  if (nrow(cache_dict[proposed_name == name, ]) == 0) {
    verify_taxon(name, itis_database, cache_dict) %>%
      data.table::fwrite(cache_path, append = TRUE)
  }
}

#' Set each unique taxon in a taxonomic dictionary a unique ID
#' 
#' Give to each unverified proposed names a different ID that ends in 0 and to 
#' each verified_name-verified_kingdom combination a different ID that ends in
#' a number corresponding to the associated kingdom.
#' 
set_taxon_ID <- function(dict, itis_database) {
  #TODO
  dict
}


# Proceed to single name verification ==========================================
#' Verify a single proposed taxonomic name (without cache)
#' 
#' Verify a single proposed taxonomic name and return its verification as well
#' as the verification of associated names in the verification process.
#' 
verify_taxon <- function(name, itis_database, cache_dict) {
  # Check the name for valid matches and synonyms in ITIS:
  itis_results_1 <- check_itis(name, itis_database)
  
  # Check the returned names for :
  gnr_results <- itis_results_1 %>%
    .[['proposed_name']] %>%
    unique() %>%
    purrr::map_df(check_gnr)
  
  # Verify the name dictionary with ITIS to fill remaining information:
  itis_results_2 <- gnr_results %>%
    .[['proposed_name']] %>%
    unique() %>%
    purrr::map_df(~check_itis(., itis_database))
  
  # Combine retrieved information:
  combined_results <- rbind(itis_results_1, itis_results_2) %>%
    merge(gnr_results, by = "proposed_name", all.x = TRUE)
  
  # Set verification information of all entries:
  combined_results %<>% set_verification_info(cache_dict)
}

#' Verify proposed name in the ITIS database
#' 
check_itis <- function(name, itis_database) {
  # Find matches in ITIS database:
  itis_matches <- itis_database$taxonomic_units[complete_name == name, ]
  
  # If entry found, process results and synonyms, else return default table:
  if (nrow(itis_matches) > 0) {
    # Left join with synonyms:
    itis_matches %<>%
      merge(itis_database$synonym_links, by = "tsn", all.x = TRUE)
    
    # Add valid synonyms:
    valid_synonyms <- itis_database$taxonomic_units[
      tsn %in% itis_matches[!is.na(tsn_accepted), ][['tsn_accepted']],
    ][
      , ':='(tsn_accepted = NA)
    ]
    
    # Bind itis_matches and valid_synonyms:
    itis_matches %<>% rbind(valid_synonyms)
    
    # Left join with kingdom and rank names:
    itis_matches %<>%
      merge(itis_database$kingdoms, by = "kingdom_id", all.x = TRUE) %>%
      merge(itis_database$ranks, by = c("kingdom_id", "rank_id"), all.x = TRUE)
    
    # Remove duplicates by choosing the entry with rank closest to species:
    itis_matches %<>% .[
      , priority := abs(220 - rank_id)
    ]
    itis_matches %<>% .[
      order(priority),
    ] %>%
      unique(by = c("complete_name", "kingdom_id"))
    
    # Format data to return:
    return_table <- itis_matches[, .(
      proposed_name = complete_name,
      found_in_itis = TRUE,
      itis_tsn = tsn,
      itis_accepted_tsn = tsn_accepted,
      itis_status = ifelse(
        n_usage %in% c("accepted", "valid"), 
        "Accepted", 
        ifelse(is.na(tsn_accepted), "Unaccepted", "Accepted synonym found")
      ),
      itis_reason = unaccept_reason,
      itis_rank = rank_name,
      itis_kingdom = kingdom_name
    )]
  } else {
    return_table <- data.table::data.table(
      proposed_name = name,
      found_in_itis = FALSE,
      itis_tsn = NA_integer_,
      itis_accepted_tsn = NA_integer_,
      itis_status = NA_character_,
      itis_reason = NA_character_,
      itis_rank = NA_character_,
      itis_kingdom = NA_character_
    )
  }
  return_table
}

#' Verify proposed name with GNR resolve
#' 
check_gnr <- function(name) {
  # Function to return default GNR result (not found):
  default_gnr <- function(name) {
    data.table::data.table(
      proposed_name = name,
      gnr_status = 'Not found',
      gnr_match = NA_character_,
      gnr_score = NA_real_,
      gnr_source = NA_character_
    )
  }
  
  # Check the name using the GNR resolve function:
  gnr_outcome <- taxize::gnr_resolve(
    name, best_match_only = "TRUE", canonical = TRUE
  )
  
  # Process the result or return default table if no result is found:
  if (nrow(gnr_outcome) > 0) {
    gnr_result <- data.table::data.table(
      proposed_name = name,
      gnr_status = ifelse(
        name == gnr_outcome$matched_name2, 
        "Perfect match",
        "Fuzzy match"
      ),
      gnr_match = gnr_outcome$matched_name2,
      gnr_score = gnr_outcome$score,
      gnr_source = gnr_outcome$data_source_title
    )
    
    # If the match is fuzzy, add a row for the fuzzy name:
    if (gnr_result$gnr_status == "Fuzzy match") {
      gnr_result %<>% rbind(data.table::data.table(
        proposed_name = gnr_outcome$matched_name2,
        gnr_status = 'Perfect match',
        gnr_match = gnr_outcome$matched_name2,
        gnr_score = 0.988,
        gnr_source = gnr_outcome$data_source_title
      ))
    }
  } else {
    gnr_result <- default_gnr(name)
  }
  gnr_result
}

#' Set verification information for all rows of a name dictionary
#' 
set_verification_info <- function(name_dict, cache_dict) {
  #TODO
  name_dict[, ':='(
    verified_name = character(),
    verified_level = character(),
    verified_rank = character(),
    verified_kingdom  = character(),
    is_verified = logical(),
    verification_status = character(),
    verification_time = character(),
    last_proposed_time = character()
  )]
  
  # Set the correct columns order
  data.table::setcolorder(name_dict, c(
    "proposed_name", "verified_name", "verified_level", "verified_rank", 
    "verified_kingdom", "is_verified", "verification_status", "gnr_status", 
    "gnr_match", "gnr_score", "gnr_source", "found_in_itis", "itis_tsn", 
    "itis_accepted_tsn", "itis_status", "itis_reason", "itis_rank", 
    "itis_kingdom", "verification_time", "last_proposed_time"
  ))
  name_dict
}
