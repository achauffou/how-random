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
    cache_dict <- data.table::fread(cache_path)
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
    lapply(check_single_name, synonyms, cache_path)
  
  # Remove duplicates, set taxon IDs and update last proposition date:
  data.table::fread(cache_path) %>%
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
    verified_kingdom  = character(),
    verified_rank = character(),
    is_verified = logical(),
    verification_status = character(),
    gnr_status = character(),
    gnr_match = character(),
    gnr_score = numeric(),
    gnr_source = character(),
    found_in_itis = character(),
    itis_tsn = integer(),
    itis_status = character(),
    itis_reason = character(),
    itis_kingdom = character(),
    itis_rank = character(),
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
  cache_dict <- data.table::fread(cache_path)
  
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
    check_gnr()
  
  # Verify the name dictionary with ITIS to fill remaining information:
  itis_results_2 <- gnr_results %>%
    .[['proposed_name']] %>%
    check_itis(itis_database)
  
  # Combine retrieved information:
  combined_results <- rbind(itis_results_1, itis_results_2) %>%
    merge(gnr_results, by = "proposed_name", all.x = TRUE)
  
  # Set verification information of all entries:
  combined_results %<>% set_verification_info(cache_dict)
}

#' Verify proposed names in the ITIS database
#' 
check_itis <- function(names, itis_database) {
  #TODO
  data.table::data.table(
    proposed_name = character(),
    found_in_itis = character(),
    itis_tsn = integer(),
    itis_status = character(),
    itis_reason = character(),
    itis_kingdom = character(),
    itis_rank = character()
  )
}

#' Verify proposed names with GNR resolve
#' 
check_gnr <- function(names) {
  #TODO
  data.table::data.table(
    proposed_name = character(),
    gnr_status = character(),
    gnr_match = character(),
    gnr_score = numeric(),
    gnr_source = character()
  )
}

#' Set verification information for all rows of a name dictionary
#' 
set_verification_info <- function(name_dict, cache_dict) {
  #TODO
  name_dict[, ':='(
    verified_name = character(),
    verified_level = character(),
    verified_kingdom  = character(),
    verified_rank = character(),
    is_verified = logical(),
    verification_status = character(),
    verification_time = character(),
    last_proposed_time = character()
  )]
}
