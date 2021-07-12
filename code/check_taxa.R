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
  cache_path
) {
  # Open cached dictionary or create one if it does not exist:
  if(file.exists(cache_path)) {
    cache_dict <- data.table::fread(
      cache_path, na.strings = c("", "NA"), colClasses = tax_dict_cols_classes()
    )
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
  if (length(names_to_verify) > 0) {
    nb_cores <- get_nb_cpus()
    message(paste("Checking", length(names_to_verify), "species names..."))
    pb <- progress::progress_bar$new(
      total = length(names_to_verify),
      format = "  :current/:total :percent |:bar| :elapsedfull",
      incomplete = " ",
      force = TRUE
    )
    new_dict <- 1:length(names_to_verify) %>%
      parallel::mclapply(
        function(x) {
          manual_kingdom <- proposed_names[
            proposed_name == names_to_verify[x], ][['manual_kingdom']]
          check_single_name(names_to_verify[x], itis_database, cache_path, 
                            manual_kingdom, nb_cores)
          if(x %% nb_cores == 0) pb$update(x/length(names_to_verify))
        }, mc.cores = nb_cores
      )
    pb$terminate()
    message()
  }

  # Remove duplicates(proposed names with same verified kingdom):
  taxonomic_dict <- data.table::fread(
    cache_path, na.strings = c("", "NA"), colClasses = tax_dict_cols_classes()
  ) %>%
    .[order(-verification_time), ] %>%
    unique(by = c("proposed_name", "verified_kingdom"))

  # Update last proposition time for names that were proposed:
  this_time <- as.character(Sys.time())
  taxonomic_dict[
    proposed_name %in% proposed_names[['proposed_name']],
    last_proposed_time := this_time
  ]

  # Save the taxonomic dictionary without duplicate lines in cache:
  data.table::fwrite(taxonomic_dict, cache_path)

  # Set genus- and species-aggregated names:
  taxonomic_dict[, ':='(
    verified_genus = stringr::str_extract(verified_name, "^(\\w*)"),
    verified_species = stringr::str_extract(verified_name, "^(\\w* ?\\w*)")
  )]
  
  # Set taxon ID to distinguish unique taxa and return the dictionary:
  taxonomic_dict %<>%
    set_taxon_ID()

  # Set columns order:
  taxonomic_dict %>%
    data.table::setcolorder(c(
      "genus_id", "sp_id", "taxon_id", "proposed_name", "verified_name", 
      "verified_genus", "verified_species", "verified_level", 
      "verified_kingdom", "verified_rank", "is_verified", "verification_status",
      "gnr_status", "gnr_match", "gnr_score", "gnr_source", "ncbi_rank", 
      "found_in_itis", "itis_tsn", "itis_accepted_tsn", "itis_status", 
      "itis_reason", "itis_rank", "found_in_pow", "pow_fqId", "pow_status", 
      "pow_accepted_name", "pow_accepted_fqId", "pow_rank", "verification_time", 
      "last_proposed_time"
    ))
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
    ncbi_rank = character(),
    found_in_itis = logical(),
    itis_tsn = integer(),
    itis_accepted_tsn = integer(),
    itis_status = character(),
    itis_reason = character(),
    itis_rank = character(),
    found_in_pow = logical(),
    pow_fqId = character(),
    pow_status = character(),
    pow_accepted_name = character(),
    pow_accepted_fqId = character(),
    pow_rank = character(),
    verification_time = character(),
    last_proposed_time = character()
  )
}

#' Get classes of the columns of a taxonomic dictionary
#'
tax_dict_cols_classes <- function() {
  list(
    character = c(1:5, 7:9, 11:12, 16:18, 20:26),
    logical = c(6, 13, 19),
    integer = c(14, 15),
    numeric = c(10)
  )
}

#' Set each unique taxon in a taxonomic dictionary a unique ID
#'
#' Give to each unverified proposed names a different ID that ends in 0 and to
#' each verified_name-verified_kingdom combination a different ID that ends in
#' a number corresponding to the associated kingdom.
#'
set_taxon_ID <- function(dict) {
  # Split dictionary into two depending on whether the name is verified
  dict_unverified <- dict[is_verified == FALSE, ]
  dict_verified <- dict[is_verified == TRUE, ]

  # Give a unique taxon ID ending in 0 to unverified species:
  dict_unverified[, ':='(taxon_id = .GRP * 10, genus_id = .GRP * 10, 
                         sp_id = .GRP * 10), by = .(proposed_name)]

  # Give a unique taxon ID to verified names:
  dict_verified[
    ,
    taxon_id := .GRP * 10 + sapply(
      verified_kingdom,
      switch,
      "Archaea" = 1,
      "Bacteria" = 2,
      "Other Eukaryota" = 3,
      "Animalia" = 4,
      "Fungus" = 5,
      "Plantae" = 6,
      7
    ),
    by = .(verified_name, verified_kingdom)
  ]
  dict_verified[
    ,
    genus_id := .GRP * 10 + sapply(
      verified_kingdom,
      switch,
      "Archaea" = 1,
      "Bacteria" = 2,
      "Other Eukaryota" = 3,
      "Animalia" = 4,
      "Fungus" = 5,
      "Plantae" = 6,
      7
    ),
    by = .(verified_genus, verified_kingdom)
  ]
  dict_verified[
    ,
    sp_id := .GRP * 10 + sapply(
      verified_kingdom,
      switch,
      "Archaea" = 1,
      "Bacteria" = 2,
      "Other Eukaryota" = 3,
      "Animalia" = 4,
      "Fungus" = 5,
      "Plantae" = 6,
      7
    ),
    by = .(verified_species, verified_kingdom)
  ]

  # Bind verified and unverified names and reorder columns:
  rbind(dict_unverified, dict_verified)
}


# Proceed to single name verification ==========================================
#' Check a single proposed taxonomic name (use and save in cache)
#'
#' Check the proposed name in the taxonomic directory in cache and verify the
#' name if it is not already present in the dictionary. If the name was not
#' already in the taxonomic dictionary, update the cache accordingly after
#' verifying it.
#'
check_single_name <- function(name, itis_database, cache_path, kingdom = NA, nb_cores = 1) {
  # Check whether the proposed name is already in the cache:
  cache_dict <- data.table::fread(
    cache_path, na.strings = c("", "NA"), colClasses = tax_dict_cols_classes()
  )

  # If the name is not already in cache, verify it:
  if (nrow(cache_dict[proposed_name == name, ]) == 0) {
    verify_taxon(name, itis_database, kingdom, nb_cores) %>%
      data.table::fwrite(cache_path, append = TRUE)
  }
}

#' Verify a single proposed taxonomic name (without cache)
#'
#' Verify a single proposed taxonomic name and return its verification as well
#' as the verification of associated names in the verification process.
#'
verify_taxon <- function(name, itis_database, kingdom = NA, nb_cores = 1) {
  # Recursively check names in several databases:
  all_names <- name
  itis_results <- NULL
  itis_names <- character()
  pow_results <- NULL
  pow_names <- character()
  gnr_results <- NULL
  gnr_names <- character()
  while (
    !all(sapply(list(itis_names, pow_names, gnr_names), identical, all_names))
  ) {
    # Check names in ITIS:
    if (length(setdiff(all_names, itis_names)) > 0) {
      itis_results <- setdiff(all_names, itis_names) %>%
        purrr::map_df(check_itis, itis_database) %>%
        rbind(itis_results)
      itis_names <- itis_results$proposed_names %>% unique() %>% sort()
    }
    
    # Check names in POW:
    if (length(setdiff(all_names, pow_names)) > 0) {
      pow_results <- setdiff(all_names, pow_names) %>%
        purrr::map_df(check_pow, nb_cores) %>%
        rbind(pow_results)
      pow_names <- pow_results$proposed_names %>% unique() %>% sort()
    }
    
    # Check names in GNR:
    if (length(setdiff(all_names, gnr_names)) > 0) {
      gnr_results <- setdiff(all_names, gnr_names) %>%
        purrr::map_df(check_gnr) %>%
        rbind(gnr_results)
      gnr_names <- gnr_results$proposed_names %>% unique() %>% sort()
    }
    
    # Update vector with all names:
    all_names <- c(itis_names, pow_names, gnr_names) %>% unique() %>% sort()
  }
  
  # Combine retrieved information:
  combined_results <- itis_results[!is.na(verified_kingdom)] %>%
    merge(pow_results[!is.na(verified_kingdom)], 
          by = c("proposed_name", "verified_kingdom"), all = TRUE) %>%
    merge(gnr_results[!is.na(verified_kingdom)], 
          by = c("proposed_name", "verified_kingdom"), all = TRUE)
  combined_results %<>% 
    merge(itis_results[is.na(verified_kingdom)], by = c("proposed_name"), all = TRUE)
  combined_results[, ':='(
    verified_kingdom = ifelse(is.na(verified_kingdom.x), verified_kingdom.y, verified_kingdom.x),
    verified_kingdom.x = NULL, verified_kingdom.y = NULL,
    found_in_itis = ifelse(is.na(found_in_itis.x), found_in_itis.y, found_in_itis.x),
    found_in_itis.x = NULL, found_in_itis.y = NULL,
    itis_tsn = ifelse(is.na(itis_tsn.x), itis_tsn.y, itis_tsn.x),
    itis_tsn.x = NULL, itis_tsn.y = NULL,
    itis_accepted_tsn = ifelse(is.na(itis_accepted_tsn.x), itis_accepted_tsn.y, itis_accepted_tsn.x),
    itis_accepted_tsn.x = NULL, itis_accepted_tsn.y = NULL,
    itis_status = ifelse(is.na(itis_status.x), itis_status.y, itis_status.x),
    itis_status.x = NULL, itis_status.y = NULL,
    itis_reason = ifelse(is.na(itis_reason.x), itis_reason.y, itis_reason.x),
    itis_reason.x = NULL, itis_reason.y = NULL,
    itis_rank = ifelse(is.na(itis_rank.x), itis_rank.y, itis_rank.x),
    itis_rank.x = NULL, itis_rank.y = NULL
  )]
  combined_results %<>% 
    merge(pow_results[is.na(verified_kingdom)], by = c("proposed_name"), all = TRUE)
  combined_results[, ':='(
    verified_kingdom = ifelse(is.na(verified_kingdom.x), verified_kingdom.y, verified_kingdom.x),
    verified_kingdom.x = NULL, verified_kingdom.y = NULL,
    found_in_pow = ifelse(is.na(found_in_pow.x), found_in_pow.y, found_in_pow.x),
    found_in_pow.x = NULL, found_in_pow.y = NULL,
    pow_fqId = ifelse(is.na(pow_fqId.x), pow_fqId.y, pow_fqId.x),
    pow_fqId.x = NULL, pow_fqId.y = NULL,
    pow_status = ifelse(is.na(pow_status.x), pow_status.y, pow_status.x),
    pow_status.x = NULL, pow_status.y = NULL,
    pow_accepted_name = ifelse(is.na(pow_accepted_name.x), pow_accepted_name.y, pow_accepted_name.x),
    pow_accepted_name.x = NULL, pow_accepted_name.y = NULL,
    pow_accepted_fqId = ifelse(is.na(pow_accepted_fqId.x), pow_accepted_fqId.y, pow_accepted_fqId.x),
    pow_accepted_fqId.x = NULL, pow_accepted_fqId.y = NULL,
    pow_rank = ifelse(is.na(pow_rank.x), pow_rank.y, pow_rank.x),
    pow_rank.x = NULL, pow_rank.y = NULL
  )]
  combined_results %<>% 
    merge(gnr_results[is.na(verified_kingdom)], by = c("proposed_name"), all = TRUE)
  combined_results[, ':='(
    verified_kingdom = ifelse(is.na(verified_kingdom.x), verified_kingdom.y, verified_kingdom.x),
    verified_kingdom.x = NULL, verified_kingdom.y = NULL,
    gnr_status = ifelse(is.na(gnr_status.x), gnr_status.y, gnr_status.x),
    gnr_status.x = NULL, gnr_status.y = NULL,
    gnr_match = ifelse(is.na(gnr_match.x), gnr_match.y, gnr_match.x),
    gnr_match.x = NULL, gnr_match.y = NULL,
    gnr_score = ifelse(is.na(gnr_score.x), gnr_score.y, gnr_score.x),
    gnr_score.x = NULL, gnr_score.y = NULL,
    gnr_source = ifelse(is.na(gnr_source.x), gnr_source.y, gnr_source.x),
    gnr_source.x = NULL, gnr_source.y = NULL,
    ncbi_rank = ifelse(is.na(ncbi_rank.x), ncbi_rank.y, ncbi_rank.x),
    ncbi_rank.x = NULL, ncbi_rank.y = NULL
  )]
  
  # Set verification information of all entries:
  combined_results %<>% set_verification_info()
  
  # Set the correct columns order:
  data.table::setcolorder(combined_results, c(
    "proposed_name", "verified_name", "verified_level", "verified_kingdom",
    "verified_rank", "is_verified", "verification_status", "gnr_status", 
    "gnr_match", "gnr_score", "gnr_source", "ncbi_rank", "found_in_itis", 
    "itis_tsn", "itis_accepted_tsn", "itis_status", "itis_reason", "itis_rank",
    "found_in_pow", "pow_fqId", "pow_status", "pow_accepted_name", 
    "pow_accepted_fqId", "pow_rank", "verification_time", "last_proposed_time"
  ))
  
  # If there is no entry for the manually specified kingdom, add one:
  if (!is.na(kingdom)) {
    if (nrow(combined_results[verified_kingdom %in% kingdom, ]) == 0) {
      combined_results <- data.table::data.table(
        proposed_name = name,
        verified_name = name,
        verified_level = c('higher', 'species', 'subspecies') %>%
          .[stringr::str_count(name, " ") + 1],
        verified_kingdom  = kingdom,
        verified_rank = NA_character_,
        is_verified = TRUE,
        verification_status = "Manual kingdom",
        gnr_status = "Not found",
        gnr_match = NA_character_,
        gnr_score = NA_real_,
        gnr_source = NA_character_,
        ncbi_rank = NA_character_,
        found_in_itis = FALSE,
        itis_tsn = NA_integer_,
        itis_accepted_tsn = NA_integer_,
        itis_status = NA_character_,
        itis_reason = NA_character_,
        itis_rank = NA_character_,
        found_in_pow = FALSE,
        pow_fqId = NA_character_,
        pow_status = NA_character_,
        pow_accepted_name = NA_character_,
        pow_accepted_fqId = NA_character_,
        pow_rank = NA_character_,
        verification_time = as.character(Sys.time()),
        last_proposed_time = NA_character_
      ) %>% rbind(combined_results)
    }
  }
  combined_results
}

#' Verify proposed name in the ITIS database
#'
check_itis <- function(name, itis_database) {
  # Find matches in ITIS database:
  itis_matches <- itis_database$taxonomic_units[complete_name == name, ]
  
  # If entry found, process results and synonyms, else return default table:
  if (nrow(itis_matches) > 0) {
    # In each kingdom, if there is a valid name, remove invalid ones:
    kingdoms_w_valid_name <- itis_matches[
      n_usage %in% c("valid", "accepted")
    ][
      , .N, by = .(kingdom_id)
    ][N > 0,][['kingdom_id']]
    itis_matches <- itis_matches[
      (n_usage %in% c("valid", "accepted")) | !(kingdom_id %in% kingdoms_w_valid_name),
    ]
    
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

    # Merge Chromista and Protozoa into one kingdom:
    itis_matches[
      kingdom_name %in% c("Protozoa", "Chromista"),
      kingdom_name := "Other Eukaryota"
    ]

    # Remove duplicates by choosing the entry with rank closest to species:
    itis_matches %<>% .[
      , priority := abs(220 - rank_id)
    ]
    itis_matches %<>% .[
      order(priority),
    ] %>%
      unique(by = c("complete_name", "kingdom_name"))

    # Format data to return:
    return_table <- itis_matches[, .(
      proposed_name = complete_name,
      verified_kingdom = kingdom_name,
      found_in_itis = TRUE,
      itis_tsn = tsn,
      itis_accepted_tsn = tsn_accepted,
      itis_status = ifelse(
        n_usage %in% c("accepted", "valid"),
        "Accepted",
        ifelse(is.na(tsn_accepted), "Unaccepted", "Synonym found")
      ),
      itis_reason = unaccept_reason,
      itis_rank = tolower(rank_name)
    )]
  } else {
    return_table <- data.table::data.table(
      proposed_name = name,
      verified_kingdom = NA_character_,
      found_in_itis = FALSE,
      itis_tsn = NA_integer_,
      itis_accepted_tsn = NA_integer_,
      itis_status = NA_character_,
      itis_reason = NA_character_,
      itis_rank = NA_character_
    )
  }
  return_table
}

#' Verify a single proposed name in Kew's Plants of the World
#' 
check_pow <- function(proposed_name, nb_cores = 1) {
  # Retrieve POW matches:
  success <- FALSE
  trial <- 0
  pow_matches <- NULL
  while(success == FALSE && trial < 10) {
    tryCatch({
      pow_matches <- taxize::pow_search(proposed_name) %$% 
        data %>% 
        data.table::as.data.table()
      success <- TRUE
    }, error = function(e) {}, warning = function(w) {}, finally = {})
    Sys.sleep(0.21 * nb_cores)
    trial <- trial + 1
  }
  
  # Return empty data table if no result is found:
  if (nrow(pow_matches) == 0) {
    return(data.table::data.table(
      proposed_name = proposed_name,
      verified_kingdom = NA_character_,
      found_in_pow = FALSE,
      pow_fqId = NA_character_,
      pow_status = NA_character_,
      pow_accepted_name = NA_character_,
      pow_accepted_fqId = NA_character_,
      pow_rank = NA_character_
    ))
  }
  
  # Format rows to be returned:
  if ("SynonymOf.accepted" %in% names(pow_matches)) {
    synonyms <- pow_matches[
      name %in% proposed_name & synonymOf.accepted == TRUE][['synonymOf.name']]
    pow_matches[, ':='(
      pow_status = ifelse(accepted == TRUE, "Accepted", ifelse(
        synonymOf.accepted == TRUE, "Synonym found", "Unaccepted")),
      pow_accepted_name = synonymOf.name,
      pow_accepted_fqId = stringr::str_remove(synonymOf.url, "/taxon/")
    )]
  } else {
    synonyms <- NULL
    pow_matches[, ':='(
      pow_status = ifelse(accepted == TRUE, "Accepted", "Unaccepted"),
      pow_accepted_name = NA_character_,
      pow_accepted_fqId = NA_character_
    )]
  }
  pow_matches[name %in% c(proposed_name, synonyms), .(
    proposed_name = proposed_name,
    verified_kingdom = kingdom,
    found_in_pow = TRUE,
    pow_fqId = fqId,
    pow_status,
    pow_accepted_name,
    pow_accepted_fqId,
    pow_rank = tolower(rank)
  )]
}

#' Verify proposed name with GNR resolve
#'
check_gnr <- function(name, nb_cores = 1) {
  # Function to return default GNR result (not found):
  default_gnr <- function(name) {
    data.table::data.table(
      proposed_name = name,
      verified_kingdom = NA_character_,
      gnr_status = 'Not found',
      gnr_match = NA_character_,
      gnr_score = NA_real_,
      gnr_source = NA_character_,
      ncbi_rank = NA_character_
    )
  }

  # Check the name using the GNR resolve function:
  gnr_outcome <- data.table::data.table(NULL)
  tryCatch({
    gnr_outcome <- taxize::gnr_resolve(
      name, best_match_only = "TRUE", canonical = TRUE
    )
  }, error = function(e) {}, warning = function(w) {}, finally = {})

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
    
    # Get the NCBI kingdom and rank:
    gnr_result %<>%
      apply(1, as.list) %>%
      lapply(get_ncbi_info, nb_cores) %>%
      data.table::rbindlist() %>%
      .[, gnr_score := as.numeric(gnr_score)]
  } else {
    gnr_result <- default_gnr(name)
  }
  gnr_result
}

#' Set verification information for all rows of a name dictionary
#'
set_verification_info <- function(name_dict) {
  # If applicable, set GNR and ITIS status to 'Not found' instead of NA:
  name_dict[is.na(found_in_itis), found_in_itis := FALSE]
  name_dict[is.na(found_in_pow), found_in_pow := FALSE]
  name_dict[is.na(gnr_status), gnr_status := "Not found"]
  
  # If there was at least one verified row, remove unverified row:
  if (nrow(name_dict[
    found_in_itis == TRUE | !gnr_status %in% "Not found" | found_in_pow == TRUE
  ]) > 0) {
    name_dict %<>% .[
      found_in_itis == TRUE | !gnr_status %in% "Not found" | found_in_pow == TRUE
    ]
  }
  
  # Set kingdom of GNR fuzzy matches without kingdom if their match has one:
  get_fuzzy_kingdom <- function(name) {
    name_dict[proposed_name %in% name][['verified_kingdom']]
  }
  name_dict <- rbind(
    name_dict[
      gnr_status %in% "Fuzzy match" & is.na(verified_kingdom), 
      .(verified_kingdom = get_fuzzy_kingdom(gnr_match)),
      by = proposed_name
    ][
      name_dict[gnr_status %in% "Fuzzy match" & is.na(verified_kingdom)][
        , verified_kingdom := NULL], 
      on = .(proposed_name)
    ],
    name_dict[!gnr_status %in% "Fuzzy match" | !is.na(verified_kingdom)]
  )
  
  # Set the same verified name for entries with identical kingdom:
  name_dict[, ':='(
    verified_name = get_verified_name(.SD)
  ), by = .(verified_kingdom)]

  # Remove abbreviations from proposed and verified names:
  name_dict[, ':='(verified_name = remove_abbreviations(verified_name))]
  dict_wo_abbr <- data.table::copy(name_dict) %>%
    .[, ':='(proposed_name = remove_abbreviations(proposed_name))]
  name_dict <- rbind(name_dict, dict_wo_abbr) %>% unique()

  # Set verified level, verification status and time for all entries:
  name_dict[, ':='(
    verified_level = c('higher', 'species', 'subspecies') %>%
      .[stringr::str_count(verified_name, " ") + 1],
    is_verified = !is.na(verified_name),
    verification_status = get_verified_status(itis_status, pow_status, gnr_status),
    verification_time = as.character(Sys.time()),
    last_proposed_time = NA_character_
  )]
  
  # Set verified rank for all entries:
  name_dict[, ':='(
    verified_rank = get_verified_rank(.SD)
  ), by = .(verified_name, verified_kingdom)]
  name_dict
}

#' Set the valid verified name of a dictionary of taxonomic synonyms
#'
get_verified_name <- function(dict) {
  # If there is an accepted ITIS name, set it as the verified name:
  valid_itis_name <- dict[itis_status == "Accepted", ][['proposed_name']][1]
  if (!is.na(valid_itis_name)) return(valid_itis_name)
  
  # Then, if there is an accepted POW name, set it as the verified name:
  valid_pow_name <- dict[pow_status == "Accepted", ][['proposed_name']][1]
  if (!is.na(valid_pow_name)) return(valid_pow_name)

  # Then, if there is a GNR match, set it as the verified name:
  valid_gnr_name <- dict[!is.na(gnr_match), ][['gnr_match']][1]
  if (!is.na(valid_gnr_name)) return(valid_gnr_name)

  # If no verified name is available, set the verified name to NA:
  return(NA_character_)
}

#' Get the verification status of a taxonomic name entry
#'
get_verified_status <- function(vec_itis_status, vec_pow_status, vec_gnr_status) {
  # Function to get status of a single entry:
  single_entry_status <- function(itis_status, pow_status, gnr_status) {
    if (!is.na(itis_status)) {
      if (itis_status == "Accepted") {
        verif_status <- "ITIS accepted"
      } else {
        verif_status <- "ITIS synonym"
      }
    } else if (pow_status %in% c("Accepted", "Synonym found")) {
      if (pow_status == "Accepted") {
        verif_status <- "POW accepted"
      } else {
        verif_status <- "POW synonym"
      }
    } else {
      verif_status <- switch(gnr_status,
                             "Not found" = "Unverified",
                             "Perfect match" = "GNR perfect match",
                             "Fuzzy match" = "GNR fuzzy match"
      )
    }
    verif_status
  }

  # Apply function to all entries:
  sapply(
    1:length(vec_itis_status),
    function(x) single_entry_status(vec_itis_status[x], vec_pow_status[x], vec_gnr_status[x])
  )
}

#' Get the verified kingdom of an entry
#' 
get_verified_rank <- function(dict) {
  # If there is an accepted ITIS rank, set it as the verified rank:
  valid_itis_rank <- dict[itis_status == "Accepted", ][['itis_rank']][1]
  if (!is.na(valid_itis_rank)) return(valid_itis_rank)
  
  # If there is an accepted POW rank, set it as the verified rank:
  valid_pow_rank <- dict[pow_status == "Accepted", ][['pow_rank']][1]
  if (!is.na(valid_pow_rank)) return(valid_pow_rank)
  
  # Then, if there is a NCBI match, set it as the verified rank:
  valid_ncbi_rank <- dict[!is.na(ncbi_rank),][
    order(-verified_level)][['ncbi_rank']][1]
  if (!is.na(valid_ncbi_rank)) return(valid_ncbi_rank)
  
  # If no verified name is available, set the verified rank to NA:
  return(NA_character_)
}

#' Get the possible NCBI kingdoms of a given dictionary row
#'
get_ncbi_info <- function(name_dict, nb_cores = 1) {
  # Convert named list to data.table:
  name_dict %<>% data.table::as.data.table()

  # Retrieve the UID of the given taxon:
  uids <- get_ncbi_uid_from_name(name_dict[['gnr_match']], nb_cores)

  # Retrieve NCBI kingdoms of the given UID:
  if (!is.null(uids)) {
    ncbi_info <- uids %>%
      lapply(get_ncbi_info_from_uid, nb_cores)
  } else {
    ncbi_info <- list(c(NA_character_, NA_character_))
  }
  
  # Return one row per possible kingdom:
  ncbi_info %>%
    purrr::map_df(~ data.table::copy(name_dict[, ':='(
      verified_kingdom = .[1], ncbi_rank = .[2]
    )]))
}

#' Get all possible UIDs of a species name by performing a NCBI eutils query
#'
get_ncbi_uid_from_name <- function(name, nb_cores = 1) {
  # Prepare query:
  api_key <- Sys.getenv("ENTREZ_KEY")
  name <- gsub(" ", "+", name)
  query <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?",
    "db=taxonomy&retmode=json&api_key=",
    api_key,
    "&term=",
    name
  )

  # Get query results:
  success <- FALSE
  trial <- 0
  sleep <- ifelse(nchar(Sys.getenv("ENTREZ_KEY")) > 0, 0.101, 0.334) * nb_cores
  uids <- NULL
  while(success == FALSE && trial < 10) {
    tryCatch({
      query_result <- query %>%
        readLines(warn = FALSE)%>%
        glue::glue_collapse() %>%
        rjson::fromJSON()
      if (length(query_result$esearchresult$ERROR) == 0) {
        uids <- query_result$esearchresult$idlist
        success <- TRUE
      }
    }, error = function(e) {}, warning = function(w) {}, finally = {})
    Sys.sleep(sleep)
    trial <- trial + 1
  }
  unlist(uids)
}

#' Return the NCBI infos (kingdom/superkingdom, rank) matching a species UID
#' 
get_ncbi_info_from_uid <- function(uid, nb_cores = 1) {
  # Retrieve NCBI kingdom/superkingdom:
  ncbi_classification <- NULL
  tryCatch({
    ncbi_classification <- taxize::classification(uid, db = "ncbi", max_tries = 10)[[1]] %>%
      data.table::as.data.table() %>%
      data.table::setnames("rank", "ncbi_rank")
    sleep <- ifelse(nchar(Sys.getenv("ENTREZ_KEY")) > 0, 0.101, 0.334) *
      (nb_cores - 1)
    Sys.sleep(sleep)
  }, error = function(e) {}, warning = function(w) {}, finally = {})
  
  if (!is.null(ncbi_classification)) {
    # Get and correct kingdom/superkingdom name to match ITIS:
    ncbi_categ <- ncbi_classification[
      ncbi_rank %in% c("kingdom", "superkingdom"), 
    ][
      order(ncbi_rank)
    ][['name']][1]
    new_categ <- ncbi_categ %>% switch(
      "Archaea" = "Archaea",
      "Bacteria" = "Bacteria",
      "Eukaryota" = "Other Eukaryota",
      "Fungi" = "Fungi",
      "Metazoa" = "Animalia",
      "Viridiplantae" = "Plantae",
      NA_character_
    )
    # Get rank:
    new_rank <- ncbi_classification[id == uid, ][['ncbi_rank']] %>%
      stringr::str_replace("^varietas$", "variety")
  } else {
    new_categ <- NA_character_
    new_rank <- NA_character_
  }
  return(c(new_categ, new_rank))
}
