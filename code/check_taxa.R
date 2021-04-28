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
  names_to_verify %>%
    sample(10) %>% # TEMPORARY!!!
    lapply(check_single_name, itis_database, cache_path)

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
  
  # Set taxon ID to distinguish unique taxa and return the dictionary:
  taxonomic_dict %>%
    set_taxon_ID()
}

#' Create an empty taxonomic dictionary
#'
empty_taxonomic_dict <- function() {
  data.table::data.table(
    proposed_name = character(),
    verified_name = character(),
    verified_level = character(),
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

#' Get classes of the columns of a taxonomic dictionary
#' 
tax_dict_cols_classes <- function() {
  list(
    character = c(1:4, 6:8, 10, 14:19),
    logical = c(5, 11),
    integer = c(12, 13),
    numeric = c(9)
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
  cols_classes <- list(
    character = c(1:5, 7:9, 11, 15:20),
    logical = c(6, 12),
    integer = c(13, 14),
    numeric = c(10)
  )
  cache_dict <- data.table::fread(
    cache_path, na.strings = c("", "NA"), colClasses = tax_dict_cols_classes()
  )

  # If the name is not already in cache, verify it:
  if (nrow(cache_dict[proposed_name == name, ]) == 0) {
    verify_taxon(name, itis_database) %>%
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
verify_taxon <- function(name, itis_database) {
  # Check the name for valid matches and synonyms in ITIS:
  itis_results_1 <- check_itis(name, itis_database)

  # Check the returned names for matches in GNR:
  gnr_results <- itis_results_1 %>%
    .[['proposed_name']] %>%
    unique() %>%
    purrr::map_df(check_gnr)
  
  # If any fuzzy match was found, repeat the ITIS check for them:
  names_to_check <- gnr_results %>%
    .[!proposed_name %in% itis_results_1[['proposed_name']], ] %>%
    .[['proposed_name']] %>%
    unique()
  if (length(names_to_check) > 0) {
    # Check names for valid matches and synonyms in ITIS:
    itis_results_2 <- names_to_check %>% 
      purrr::map_df(~check_itis(., itis_database))
    
    # If any synonym was found, repeat the GNR check for them:
    names_to_check <- itis_results_2 %>%
      .[!proposed_name %in% gnr_results[['proposed_name']], ] %>%
      .[['proposed_name']] %>%
      unique()
    if (length(names_to_check) > 0) {
      gnr_results <- names_to_check %>% 
        purrr::map_df(check_gnr) %>% 
        rbind(gnr_results)
    }
  } else {
    itis_results_2 <- data.table::data.table(
      proposed_name = character(),
      found_in_itis = character(),
      itis_tsn = integer(),
      itis_accepted_tsn = integer(),
      itis_status = character(),
      itis_reason = character(),
      itis_rank = character(),
      itis_kingdom = character()
    )
  }
  
  # Combine retrieved information:
  combined_results <- rbind(itis_results_1, itis_results_2) %>%
    merge(gnr_results, by = "proposed_name", all.x = TRUE)
  
  # Set verification information of all entries:
  combined_results %<>% set_verification_info()
  
  # Set the correct columns order
  data.table::setcolorder(combined_results, c(
    "proposed_name", "verified_name", "verified_level", "verified_kingdom", 
    "is_verified", "verification_status", "gnr_status", "gnr_match", 
    "gnr_score", "gnr_source", "found_in_itis", "itis_tsn", 
    "itis_accepted_tsn", "itis_status", "itis_reason", "itis_rank",
    "itis_kingdom", "verification_time", "last_proposed_time"
  ))
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
      found_in_itis = TRUE,
      itis_tsn = tsn,
      itis_accepted_tsn = tsn_accepted,
      itis_status = ifelse(
        n_usage %in% c("accepted", "valid"),
        "Accepted",
        ifelse(is.na(tsn_accepted), "Unaccepted", "Synonym found")
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
set_verification_info <- function(name_dict) {
  # Verify the kingdom for all entries with GNR match found in ITIS:
  name_dict[
    gnr_status != "Not found", 
    verified_kingdom := .SD[order(itis_kingdom), ][['itis_kingdom']][1], 
    by = .(gnr_match)
  ]
  
  # Verify the kingdom for entries found only in GNR:
  name_dict <- rbind(
    name_dict[gnr_status == "Not found" || !is.na(verified_kingdom), ],
    name_dict[gnr_status != "Not found" && is.na(verified_kingdom), ] %>%
      apply(1, as.list) %>%
      lapply(get_ncbi_kingdoms) %>%
      data.table::rbindlist()
  )
  
  # Set the same verified name for entries with identical kingdom:
  name_dict[
    , ':='(verified_name = get_verified_name(.SD)), by = .(verified_kingdom)
  ]
  
  # Set verified level, verification status and time for all entries:
  name_dict[, ':='(
    verified_level = c('higher', 'species', 'subspecies') %>%
      .[stringr::str_count(verified_name, " ") + 1],
    is_verified = !is.na(verified_name),
    verification_status = get_verified_status(itis_status, gnr_status),
    verification_time = as.character(Sys.time()),
    last_proposed_time = NA_character_
  )]
  name_dict
}

#' Set the valid verified name of a dictionary of taxonomics synonyms
#' 
get_verified_name <- function(dict) {
  # If there is an accepted ITIS name, set it as the verified name:
  valid_itis_name <- dict[itis_status == "Accepted", ][['proposed_name']][1]
  if (!is.na(valid_itis_name)) return(valid_itis_name)
  
  # Then, if there is a GNR match, set it as the verified name:
  valid_gnr_name <- dict[!is.na(gnr_match), ][['gnr_match']][1]
  if (!is.na(valid_gnr_name)) return(valid_gnr_name)
  
  # If no verified name is available, set the verified name to NA:
  return(NA)
}

#' Get the verification status of a taxonomic name entry
#' 
get_verified_status <- function(vec_itis_status, vec_gnr_status) {
  # Function to get status of a single entry:
  single_entry_status <- function(itis_status, gnr_status) {
    if (!is.na(itis_status)) {
      if (itis_status == "Accepted") {
        verif_status <- "ITIS accepted"
      } else {
        verif_status <- "ITIS synonym"
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
    function(x) single_entry_status(vec_itis_status[x], vec_gnr_status[x])
  )
}

#' Get the possible NCBI kingdoms of a given dictionary row 
#' 
get_ncbi_kingdoms <- function(name_dict) {
  # Convert named list to data.table:
  name_dict %<>% data.table::as.data.table()
  
  # Retrieve the UID of the given taxon:
  uids <- get_ncbi_uid_from_name(name_dict[['gnr_match']])
  
  # Retrieve NCBI kingdoms of the given UID:
  if (!is.null(uids)) {
    kingdoms <- uids %>%
      sapply(get_ncbi_kingdom_from_uid)
  } else {
    kingdoms <- "Unknown"
  }
  
  # Return one row per possible kingdom:
  kingdoms %>% 
    purrr::map_df(~ data.table::copy(name_dict[, verified_kingdom := .]))
}

#' Get all possible UIDs of a species name by performing a NCBI eutils query
#' 
get_ncbi_uid_from_name <- function(name) {
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
  sleep <- ifelse(nchar(Sys.getenv("ENTREZ_KEY")) > 0, 0.101, 0.334)
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

#' Return all NCBI kingdom (or superkingdom) matching a species UID
#' 
get_ncbi_kingdom_from_uid <- function(uid) {
  # Retrieve NCBI kingdom/superkingdom:
  ncbi_categ <- taxize::classification(uid, db = "ncbi", max_tries = 10)[[1]] %>%
    data.table::as.data.table() %>%
    data.table::setnames("rank", "ncbi_rank") %>%
    .[ncbi_rank %in% c("kingdom", "superkingdom"), ] %>%
    .[order(ncbi_rank)]
  ncbi_categ <- ncbi_categ[['name']][1]
  
  # Correct kingdom/superkingdom name to match ITIS:
  if (!is.null(ncbi_categ)) {
    new_categ <- ncbi_categ %>% switch(
      "Archaea" = "Archaea",
      "Bacteria" = "Bacteria",
      "Eukaryota" = "Other Eukaryota",
      "Fungi" = "Fungi",
      "Metazoa" = "Animalia",
      "Viridiplantae" = "Plantae",
      "Unknown"
    )
  } else {
    new_categ <- "Unknown"
  }
  new_categ
}
