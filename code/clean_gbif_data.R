# Process GBIF raw archives if necessary =======================================
#' Extract and clean only GBIF archives that have not been processed yet
#' 
process_gbif_raw_archives <- function(
  archive_paths, dest_folder, land_data, gbif_keys, cache_path, 
  db_file = "cleaned_data.sqlite", table_name = "occurrences", stats_name = "cleaning_stats.csv", 
  chunk_size = getOption("CHUNK_SIZE", default = 5E4)
) {
  # If the destination folder does not exist yet, create it:
  db_path <- file.path(dest_folder, db_file)
  dir.create(dirname(db_path), showWarnings = FALSE, recursive = TRUE)
  
  # Read GBIF processing cache:
  cache <- update_gbif_processing_cache(cache_path)
  
  # Open connection with database:
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = db_path)
  on.exit(RSQLite::dbDisconnect(db))
  db_modified <- FALSE
  
  # Create the SQL database if it does not exist yet:
  if (!RSQLite::dbExistsTable(db, table_name)) {
    RSQLite::dbCreateTable(db, table_name, get_gbif_clean_occurrences_fields())
    RSQLite::dbWriteTable(
      db, "last_occurrences_modif", data.frame(time = as.character(Sys.time()))
    )
  }
  
  # Determine archives that must be extracted and processed:
  archives_processed <- cache[
    dest_folder %in% dest_folder & !is.na(processing_time)
  ][['archive_path']]
  archives_processed_names <- archives_processed %>%
    basename() %>%
    tools::file_path_sans_ext()
  archives_to_process <- setdiff(archive_paths, archives_processed)
  archives_to_extract <- setdiff(
    archive_paths, cache[dest_folder %in% dest_folder][['archive_path']]
  )
  
  # Create cleaning statistics csv or remove its irrelevant rows:
  stats_path <- file.path(dest_folder, stats_name)
  if (!file.exists(stats_path)) {
    data.table::fwrite(empty_cleaning_stats(), stats_path)
  } else {
    data.table::fread(
      stats_path, na.strings = c("", "NA"), 
      colClasses = list(character = 1:4)
    ) %>% 
      .[archive_name %in% archives_processed_names] %>% 
      data.table::fwrite(stats_path)
  }
  
  # Remove cleaned occurrences from outdated datasets if there are:
  archives_in_db <- RSQLite::dbGetQuery(db, paste(
    "SELECT DISTINCT archive_name FROM", table_name
  ))[['archive_name']]
  if (!all(archives_in_db %in% archives_processed_names)) {
    RSQLite::dbExecute(db, paste(
      "DELETE FROM", table_name, "WHERE archive_name NOT IN (", paste(
        archives_processed_names, collapse = ","
      ), ")"
    ))
    RSQLite::dbRemoveTable(db, "last_occurrences_modif")
    RSQLite::dbWriteTable(
      db, "last_occurrences_modif", data.frame(time = as.character(Sys.time()))
    )
  }
  
  # Extract archives that need to be extracted:
  if (length(archives_to_extract) > 0) {
    nb_cores <- parallel::detectCores()
    total_size <- sum(file.info(archives_to_extract)$size) / 1073741824
    message(paste("Extracting", format(total_size, digits = 4), 
                  "GB of raw GBIF archives..."))
    archives_to_extract %>%
      parallel::mclapply(function (x) {
        extract_single_gbif_archive(x, dest_folder, cache_path)
      }, mc.cores = nb_cores)
  }
  
  # Clean and save to database archives that need to be processed:
  if (length(archives_to_process) > 0) {
    # Clean and process archives one by one:
    for (archive_path in c(archives_to_process,archives_to_extract)) {
      extracted_path <- file.path(
        dest_folder, paste0(basename(tools::file_path_sans_ext(archive_path)), ".txt")
      )
      clean_gbif_occurrences(
        extracted_path, db, table_name, land_data, gbif_keys, stats_path, 
        chunk_size
      )
      cache <- data.table::fread(
        cache_path, na.strings = c("", "NA"), 
        colClasses = list(character = 1:4)
      )
      cache[archive_path %in% archive_path & dest_folder %in% dest_folder, ':='(
        processing_time = as.character(Sys.time())
      )]
      data.table::fwrite(cache, cache_path)
    }
    
    # Remove duplicate rows from database:
    RSQLite::dbExecute(db, paste(
      "DELETE FROM", table_name, "WHERE rowid NOT IN (SELECT MIN(rowid) FROM",
      table_name, 
      "GROUP BY archive_name, verified_id, decimalLatitude, decimalLongitude)"
    ))
    
    # Update last modification time of database:
    RSQLite::dbRemoveTable(db, "last_occurrences_modif")
    RSQLite::dbWriteTable(
      db, "last_occurrences_modif", data.frame(time = as.character(Sys.time()))
    )
  }
  
  # Return the database last modification date:
  RSQLite::dbReadTable(db, "last_occurrences_modif")[['time']]
}

#' Return an empty GBIF processing cache
#' 
empty_gbif_processing_cache <- function() {
  data.table::data.table(
    archive_path = character(),
    dest_folder = character(),
    extraction_time = character(),
    processing_time = character()
  )
}

#' Update the GBIF processing cache
#' 
update_gbif_processing_cache <- function(cache_path) {
  # Open cache file for GBIF processing or create it if it does not exist:
  if(file.exists(cache_path)) {
    cache <- data.table::fread(
      cache_path, na.strings = c("", "NA"), 
      colClasses = list(character = 1:4)
    )
  } else {
    dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
    cache <- empty_gbif_processing_cache()
  }
  
  # Remove any duplicate:
  cache %<>% unique(by = c("archive_path", "dest_folder"))
  
  # Remove processed datasets that are older than their archive:
  is_outdated <- function(archive_path, extraction_time) {
    archive_time <- as.POSIXct(file.info(archive_path)$mtime)
    archive_time > extraction_time
  }
  cache <- cache[!purrr::map2_lgl(archive_path, extraction_time, is_outdated)]
  
  # If it was not extracted, check that extraction file still exists:
  cache <- cache[!is.na(processing_time) | file.exists(file.path(
    dest_folder, 
    paste0(basename(tools::file_path_sans_ext(archive_path)), ".txt")
  ))]
  
  # Save updated cache:
  data.table::fwrite(cache, cache_path)
  cache
}


# Extract raw GBIF occurrences =================================================
#' Extract a single GBIF archive and update the raw and update the cache
#' 
extract_single_gbif_archive <- function(archive_path, dest_folder, cache_path) {
  # Destination file:
  dest_file <- file.path(
    dest_folder, paste0(basename(tools::file_path_sans_ext(archive_path)), ".txt")
  )
  
  # Extract archive to destination file:
  system(paste("unzip -pq", archive_path, "occurrence.txt > ", dest_file))
  
  # Update cache:
  data.table::data.table(
    archive_path = archive_path,
    dest_folder = dest_folder,
    extraction_time = as.character(Sys.time()),
    processing_time = NA_character_
  ) %>% data.table::fwrite(cache_path, append = TRUE)
}

#' Extract several raw GBIF archives to a combined csv file
#'
extract_gbif_archives <- function(
  raw_archives, dest_file, verbose = TRUE
) {
  # Create directory of destination file if it does not exist:
  dir.create(dirname(dest_file), showWarnings = FALSE, recursive = TRUE)
  
  # Extract first archive to the destination:
  message("Extracting GBIF occurrence data... Can take up to several hours!")
  system(paste("unzip -pq", raw_archives[1], "occurrence.txt > ", dest_file))
  
  # Extract rows from other archives and append them to the file:
  for (archive in raw_archives[-1]) {
    system(paste(
      "unzip -qp", raw_archives[1], "occurrence.txt | tail -n+2 >> ", dest_file
    ))
  }
  Sys.time()
}


# Clean GBIF occurrences =======================================================
#' Clean raw GBIF occurrences chunk by chunk and store them in a database
#' 
clean_gbif_occurrences <- function(
  occurrences_file, db, table_name, land_data, gbif_keys, stats_path, 
  chunk_size = 1E5
) {
  # Clean and assign to organisms all occurrences chunk by chunk:
  archive_name <- occurrences_file %>%
    basename() %>%
    tools::file_path_sans_ext()
  file_size <- file.info(occurrences_file) %>% .[['size']] / 1073741824
  message(paste("Cleaning", format(file_size, digits = 4), "GB of GBIF occurrences..."))
  suppressWarnings({
    suppressMessages({
      readr::read_delim_chunked(
        occurrences_file,
        function(x, idx) clean_gbif_occurrences_chunk(
          x, db, table_name, land_data, gbif_keys, stats_path, archive_name
        ),
        chunk_size = chunk_size,
        delim = "\t",
        quote = "",
        trim_ws = TRUE,
        progress = TRUE
      )
    })
  })
}

#' Return the fields of the database table that stores GBIF occurrences
#' 
get_gbif_clean_occurrences_fields <- function() {
  list(
    archive_name="TEXT",
    verified_id="INTEGER",
    decimalLatitude="REAL",
    decimalLongitude="REAL"
  )
}

#' Return empty cleaning stats file
#' 
empty_cleaning_stats <- function() {
  data.table::data.table(
    archive_name = character(),
    initial_nb = integer(),
    nb_wo_crucial_na = integer(),
    nb_wo_zero_individuals = integer(),
    nb_w_acceptable_uncertainty = integer(),
    nb_coord_invalid = integer(),
    nb_coord_capitals = integer(),
    nb_coord_centroid = integer(),
    nb_coord_equal = integer(),
    nb_coord_gbif = integer(),
    nb_coord_inst = integer(),
    nb_coord_sea = integer(),
    nb_w_clean_coordinates = integer(),
    nb_matching_gbif_key = integer(),
    nb_unique_rows = integer(),
    cleaning_time_in_secs = numeric()
  )
}


# Clean a chunk of GBIF occurrences ============================================
#' Clean and append to database a single chunk of GBIF raw data
#' 
clean_gbif_occurrences_chunk <- function(
  chunk, db, table_name, land_data, gbif_keys, stats_path, archive_name
) {
  # Start cleaning time:
  start_time <- Sys.time()
  
  # Convert chunk to data.table:
  chunk %<>% data.table::as.data.table()
  initial_nb = nrow(chunk)
  
  # Format and keep only relevant columns:
  chunk <- chunk[, .(
    taxonKey = as.integer(as.character(taxonKey)),
    genusKey = as.integer(as.character(genusKey)),
    speciesKey = as.integer(as.character(speciesKey)),
    decimalLatitude = as.numeric(as.character(decimalLatitude)),
    decimalLongitude = as.numeric(as.character(decimalLongitude)),
    coordUncert = as.numeric(as.character(coordinateUncertaintyInMeters)),
    countryCode = countrycode::countrycode(countryCode, origin = "iso2c", 
                                           destination = "iso3c"),
    year = as.integer(as.character(year)),
    indivCount = as.integer(individualCount),
    basisOfRecord,
    datasetKey,
    taxonRank
  )]
  
  # Remove rows with NAs in crucial fields:
  chunk <- na.omit(chunk, cols = c("decimalLongitude", "decimalLatitude", 
                                   "countryCode", "indivCount"))
  nb_wo_crucial_na <- nrow(chunk)
  
  # Remove rows that have zero as individuals count:
  chunk <- chunk[indivCount > 0 | is.na(indivCount)]
  nb_wo_zero_individuals <- nrow(chunk)
  
  # Remove rows that have too high coordinates uncertainty:
  chunk <- chunk[coordUncert/1000 <= 100 | is.na(coordUncert)]
  nb_w_acceptable_uncertainty <- nrow(chunk)
  
  # Clean coordinates:
  if (nb_w_acceptable_uncertainty > 0) {
    suppressWarnings({chunk <- chunk %>% 
      CoordinateCleaner::clean_coordinates(
        tests = c("capitals", "centroids", "equal", "gbif", "institutions", 
                  "zeros", "seas"),
        lon = "decimalLongitude",
        lat = "decimalLatitude",
        species = "taxonKey",
        countries = "countryCode",
        seas_ref = land_data,
        seas_scale = 10,
        verbose = FALSE
      ) %>% data.table::as.data.table()})
    nb_coord_capitals <- nrow(chunk[.cap == FALSE])
    nb_coord_centroid <- nrow(chunk[.cen == FALSE])
    nb_coord_equal <- nrow(chunk[.equ == FALSE])
    nb_coord_gbif <- nrow(chunk[.gbf == FALSE])
    nb_coord_inst <- nrow(chunk[.inst == FALSE])
    nb_coord_sea <- nrow(chunk[.sea == FALSE])
    nb_coord_invalid <- nrow(chunk[.val == FALSE])
    chunk <- chunk[.summary == TRUE]
    nb_w_clean_coordinates <- nrow(chunk)
    
    # Attribute rows to their organisms:
    chunk %<>%
      data.table::melt(
        id.vars = c("decimalLongitude", "decimalLatitude"),
        measure.vars = c("taxonKey", "genusKey", "speciesKey"),
        value.name = "gbif_key",
        na.rm = TRUE
      ) %>% 
      .[, .(gbif_key, decimalLongitude, decimalLatitude)]
    nb_matching_gbif_key <- nrow(chunk)
    chunk %<>%
      unique() %>%
      merge(gbif_keys[, .(verified_id, gbif_key)], by = "gbif_key") %>%
      .[, .(archive_name = archive_name, verified_id, decimalLatitude, decimalLongitude)]
    nb_unique_rows <- nrow(chunk)
    
    # Append cleaned chunk to the occurrences database:
    RSQLite::dbAppendTable(db, table_name, chunk)
  } else {
    nb_coord_capitals <- 0
    nb_coord_centroid <- 0
    nb_coord_equal <- 0
    nb_coord_gbif <- 0
    nb_coord_inst <- 0
    nb_coord_sea <- 0
    nb_coord_invalid <- 0
    nb_w_clean_coordinates <- 0
    nb_matching_gbif_key <- 0
    nb_unique_rows <- 0
  }
  
  # Append cleaning stats for this chunk to file:
  end_time <- Sys.time()
  time_diff <- difftime(end_time, start_time, units == "secs")
  data.table::data.table(
    archive_name = archive_name,
    initial_nb = initial_nb,
    nb_wo_crucial_na = nb_wo_crucial_na,
    nb_wo_zero_individuals = nb_wo_zero_individuals,
    nb_w_acceptable_uncertainty = nb_w_acceptable_uncertainty,
    nb_coord_invalid = nb_coord_invalid,
    nb_coord_capitals = nb_coord_capitals,
    nb_coord_centroid = nb_coord_centroid,
    nb_coord_equal = nb_coord_equal,
    nb_coord_gbif = nb_coord_gbif,
    nb_coord_inst = nb_coord_inst,
    nb_coord_sea = nb_coord_sea,
    nb_w_clean_coordinates = nb_w_clean_coordinates,
    nb_matching_gbif_key = nb_matching_gbif_key,
    nb_unique_rows = nb_unique_rows,
    cleaning_time_in_secs = as.numeric(time_diff)
  ) %>% data.table::fwrite(stats_path, append = TRUE)
}
