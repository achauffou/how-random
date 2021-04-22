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
  synonyms_db, 
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
  
  # Return the updated taxonomic dictionary with new last proposition date:
  data.table::fread(cache_path) %>%
    merge_duplicate_propositions(proposed_names[['proposed_name']]) %>%
    .[proposed_name %in% proposed_names, last_proposed_time := character(Sys.time())]
}

#' Create an empty taxonomic dictionary
#' 
empty_taxonomic_dict <- function() {
  data.table::data.table(
    taxon_ID = integer(),
    proposed_name = character(),
    proposed_level = character(),
    verified_tsn = integer(),
    verified_name = character(),
    verified_level = character(),
    verified_kingdom  = character(),
    is_verified = logical(),
    verification_status = character(),
    found_in_itis = logical(),
    itis_status = character(),
    itis_reason = character(),
    gnr_match = character(),
    gnr_score = numeric(),
    gnr_source = character(),
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
check_single_name <- function(name, synonyms_db, cache_path) {
  # Check whether the proposed name is already in the cache:
  cache_dict <- data.table::fread(cache_path)
  
  # If the name is not already in cache, verify it:
  if (nrow(cache_dict[proposed_name == name, ]) == 0) {
    verify_taxon(name, synonyms_db, cache_dict) %>%
      data.table::fwrite(cache_path, append = TRUE)
  }
}

#' Safely merge duplicate proposed names from a taxonomic dictionary
#' 
#' If there are duplicate entries for a proposed name in the given taxonomic 
#' dictionary, this function safely merge them by ensuring that their related 
#' names (that shared taxon_id with the proposed name) are given a new name. In
#' case a name from the vector of proposed names is absent, the function yields
#' an error to prevent missing verifications.
#' 
merge_duplicate_propositions <- function(dict, proposed_names) {
  #TODO
  dict
}

# Proceed to single name verification ==========================================
#' Verify a single proposed taxonomic name (without cache)
#' 
#' Verify a single proposed taxonomic name and return its verification as well
#' as the verification of associated names in the verification process.
#' 
verify_taxon <- function(name, synonyms_db, cache_dict) {
  #TODO
  empty_taxonomic_dict()
}
