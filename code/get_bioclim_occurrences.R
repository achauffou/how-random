# Extract and stack bioclimatic rasters ========================================
#' Extract and stack together raster files from multiple archives
#' 
stack_bioclim_archives <- function(
  archives, extent_coords, save_path, 
  temp_dir = file.path(tempdir(), "bioclim"), download_date = NA
) {
  # Create output folder if it does not exist:
  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
  
  # Skip the task if the rasters have already been extracted and are up-to-date:
  skip <- FALSE
  raster_files <- list.files(save_path, full.names = T)
  if (length(raster_files) > 0) {
    if (max(file.info(archives)$ctime) < min(file.info(raster_files)$mtime)) {
      skip <- TRUE
    }
  }
  
  # Extract archives to temporary directory and crop them:
  if (!skip) {
    lapply(archives, unzip, exdir = temp_dir)
    suppressWarnings({list.files(temp_dir, full.names = T) %>% 
      lapply(crop_raster_from_file, extent_coords, save_path)})
  }
  
  # Stack together all rasters:
  raster::stack(list.files(save_path, full.names = T))
}

#' Read a raster file, crop it to a given extent and save it
#' 
crop_raster_from_file <- function(file_path, extent_coords, save_path) {
  out <- raster::stack(file_path)
  raster::crs(out) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  ext <- as(raster::extent(extent_coords), 'SpatialPolygons')
  raster::crs(ext) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  ext <- raster::crop(out, ext)
  layer_name <- names(ext) %>% stringr::str_replace("^.*(2\\.5m\\_|arcmin\\_)", "")
  names(ext) <- layer_name
  file_path <- file.path(save_path, paste0(layer_name, ".tif"))
  raster::writeRaster(ext, filename = file_path, format = "GTiff", overwrite = TRUE)
  file_path
}


# Get connection to bioclimatic cache database =================================
#' Get connection of the cache database for bioclimatic conditions
#' 
#' This function returns a database connection to the cache database for 
#' bioclimatic conditions. It checks that the database for given buffers is 
#' up-to-date with stack files and (re)creates it if necessary.
#' 
get_bioclim_cache_db_con <- function(
  cache_folder, stack_path, brick, buffers = c(5000, 10000)
) {
  # Order buffers:
  buffers %<>% unique() %>% sort()
  
  # Create cache database if it does not exist or is outdated:
  db_path <- paste0(c(0, buffers), collapse = "_") %>%
    paste("bioclim_vars_cache", ., sep = "_") %>%
    paste0(".sqlite") %>%
    file.path(cache_folder, .)
  if (file.exists(db_path)) {
    if (file.info(db_path)$ctime > file.info(stack_path)$ctime) {
      db <- RSQLite::dbConnect(RSQLite::SQLite(), db_path)
    } else {
      file.remove(db_path)
      db <- create_bioclim_cache_db(db_path, brick)
    }
  } else {
    db <- create_bioclim_cache_db(db_path, brick)
  }
  db
}

#' Create a cache database for bioclimatic variables
#' 
create_bioclim_cache_db <- function(
  db_path, brick, table_name = "bioclim_vars"
) {
  # Create cache folder if it does not exist:
  dir.create(dirname(db_path), showWarnings = FALSE, recursive = TRUE)
  
  # Create the database:
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = db_path)
  paste0(names(brick), " REAL", collapse = ", ") %>%
    paste0("CREATE TABLE ", table_name, 
           " (cell INTEGER PRIMARY KEY, buffer REAL, ", ., ")") %>%
    RSQLite::dbExecute(db, .)
  db
}


# Retrieve thinned bioclimatic conditions at given coordinates =================
#' Return thinned bioclimatic conditions for a set of coordinates (cache)
#' 
#' Given a list with occurrences coordinates, this function returns a data.table
#' with the bioclimatic conditions at those coordinates. The function first 
#' thins the set of coordinates to return only one entry per grid cell of the 
#' brick. Then, it for each cell included, the function retrieves the 
#' bioclimatic conditions from the raster brick at the given location (using a 
#' cache to ensure that bioclimatic conditions of any cell are not retrieved 
#' several times). See the details of how bioclimatic conditions are retrieved 
#' at a location in the extract_cell_bioclim() function. An additional buffer 
#' column is returned with the buffer radius used to retrieve conditions.
#' 
get_thinned_bioclim_w_cache <- function(
  occurrences, brick, db, table_name = "bioclim_vars", 
  buffers = c(5000, 10000)
) {
  # Get thinned occurrences and retrieve bioclimatic conditions:
  thin_occurrences(occurrences, brick) %$%
    cell %>%
    lapply(extract_cell_bioclim, brick, brick, db, buffers) %>%
    data.table::rbindlist()
}

#' Return thinned bioclimatic conditions for a set of coordinates (no cache)
#' 
#' Given a list with occurrences coordinates, this function returns a data.table
#' with the bioclimatic conditions at those coordinates. The function first 
#' thins the set of coordinates to return only one entry per grid cell of the 
#' brick. Then, it for each cell included, the function retrieves the 
#' bioclimatic conditions from the raster brick at the given location (without 
#' using cache information). See the details of how bioclimatic conditions are 
#' retrieved at a location in the extract_cell_bioclim() function. An 
#' additional buffer column is returned with the buffer radius used to retrieve 
#' conditions. The function retrieves in a single call bioclimatic variables 
#' for all cells, which is much faster than its cached counterpart.
#' 
get_thinned_bioclim_wo_cache <- function(
  occurrences, brick, buffers = c(5000, 10000)
) {
  # Get thinned occurrences and retrieve bioclimatic conditions:
  thin_occurrences(occurrences, brick) %$%
    cell %>%
    extract_cells_bioclim(brick, buffers)
}

#' Get the unique brick cells within which occurrences fall
#' 
thin_occurrences <- function(occurrences, brick) {
  thinned <- unique(raster::cellFromXY(brick, occurrences)) %>%
    data.table::data.table(cell = .)
  thinned[!is.na(cell)]
  thinned
}


# Get bioclimatic conditions of a cell =========================================
#' Get bioclimatic conditions for a given cell (using cache database)
#' 
get_cell_bioclim <- function(cell, brick, db, table_name, buffers) {
  # # Attempt to retrieve the entry directly from the cache:
  entry <- data.table::data.table()
  tryCatch({
    entry <- paste0("SELECT * FROM ", table_name, " WHERE cell=", cell) %>%
      RSQLite::dbGetQuery(db, ., "PRAGMA busy_timeout = 10 * 1000") %>%
      data.table::as.data.table()
  }, error = function(e) {}, finally = {})

  # If an entry was found in the cache, return it:
  if (nrow(entry) > 0) {
    return(entry[1,])
  }
  
  # If the entry does not exist in cache, add it and return it:
  entry <- extract_cell_bioclim(cell, brick, buffers)
  tryCatch({
    RSQLite::dbAppendTable(db, table_name, entry, "PRAGMA busy_timeout = 10 * 1000")
  }, error = function(e) {}, finally = {})
  entry
}

#' Extract bioclimatic conditions of a given cell from a brick
#' 
#' If for any bioclimatic condition there is no value (NA) at a querried 
#' location, the average value in given buffers of increasingly larger radius 
#' are retrieved. Failing to retrieve a value (because all surrounding cells in 
#' the largest buffer have no value), the bioclimatic variable is set to NA and
#' the buffer is set to -1.
#' 
extract_cell_bioclim <- function(cell, brick, buffers) {
  # Try retrieving values at the exact location:
  buffer <- 0.0
  vars <- raster::extract(brick, cell)
  
  # If there were NAs in the variables, try using buffers
  if (any(is.na(vars))) {
    var_names <- names(vars)
    for (buffer in buffers) {
      new_vars <- raster::xyFromCell(brick, cell) %>%
        raster::extract(
          brick[[which(c(is.na(vars)))]], ., buffer = buffer, 
          cellnumbers = TRUE, fun = function(x) mean(x, na.rm = TRUE)
        )
      vars[is.na(vars)] <- new_vars
      if (all(!is.na(vars))) break
    }
    if (any(is.na(vars))) buffer = -1.0
  }
  
  # Create and return a data.table with the variables:
  data.table::data.table(cell = cell, buffer = buffer) %>%
    cbind(data.table::as.data.table(vars))
}

#' Extract bioclimatic conditions of several cells from a brick
#' 
extract_cells_bioclim <- function(cells, brick, buffers) {
  # Try retrieving values at the exact location:
  vars <- raster::extract(brick, cells)
  cells_buffer <- rep(0.0, times = length(cells))
  
  # If there were NAs in the variables, try using buffers:
  na_rows <- which(apply(vars, 1, function(x) {any(is.na(x))}))
  if (length(na_rows) > 0) {
    na_vars <- which(is.na(vars[na_rows, ]))
    for (buffer in buffers) {
      new_vars <- raster::xyFromCell(brick, cells[na_rows]) %>%
        raster::extract(
          brick, ., buffer = buffer, 
          cellnumbers = TRUE, fun = function(x) mean(x, na.rm = TRUE)
        )
      cells_buffer[na_rows] <- buffer
      vars[na_rows, ][na_vars] <- new_vars[na_vars]
      na_rows <- which(apply(vars, 1, function(x) {any(is.na(x))}))
      if (length(na_rows) == 0) break
      na_vars <- which(is.na(vars[na_rows, ]))
    }
    cells_buffer[na_rows] <- -1
  }
  
  # Create and return a data.table with the variables:
  data.table::data.table(cell = cells, buffer = cells_buffer) %>%
    cbind(data.table::as.data.table(vars))
}


# Thin and retrieve bioclimatic conditions of all GBIF occurrences =============
#' Get GBIF entities to thin
#' 
get_gbif_entities_to_thin <- function(
  gbif_keys, db_folder, db_file = "occurrences.sqlite", table_name = "cleaned",
  last_occ_update = NA
) {
  # Open connection to the database:
  db_path <- file.path(db_folder, db_file)
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = db_path)
  on.exit(RSQLite::dbDisconnect(db))
  
  # Prepare query to retrieve the entities to thin:
  keys <- gbif_keys[['gbif_key']] %>% 
    unique() %>% 
    sort() %>% 
    paste0(collapse = ", ") %>%
    paste0("(", ., ")", collapse = )
  query <- paste0(
      "SELECT DISTINCT genusKey, speciesKey, taxonKey FROM ", table_name,
      " WHERE genusKey IN ", keys, " OR speciesKey IN ", keys, 
      " OR taxonKey IN ", keys, " ORDER BY genusKey, speciesKey, taxonKey"
    )
  
  # Execute query and return outcome as a data.table:
  RSQLite::dbGetQuery(db, query) %>%
    data.table::as.data.table()
}

#' Thin and retrieve bioclimatic variables for all given GBIF entities (cached)
#' 
thin_retrieve_gbif_entities <- function(
  entities, occ_folder, stack, cache_folder, stack_path, 
  thinned_db_file = "bioclim.sqlite", occ_file = "occurrences.sqlite", 
  occ_table = "cleaned", use_raster_brick = TRUE, buffers = c(5000, 10000)
) {
  # Order buffers:
  buffers %<>% unique() %>% sort()
  
  # Check if the entities database up-to-date and create it if necessary:
  thinned_path <- file.path(occ_folder, thinned_db_file)
  thinned_db <- get_gbif_thinned_db_con(thinned_path, stack_path, stack)
  on.exit(RSQLite::dbDisconnect(thinned_db))
  
  # Open connection with GBIF occurrences database:
  occ_path <- file.path(occ_folder, occ_file)
  occ_db <- RSQLite::dbConnect(RSQLite::SQLite(), occ_path)
  on.exit(RSQLite::dbDisconnect(occ_db), add = TRUE)
  
  # Create entities cache if it does not exist or is outdated:
  thinned_cache_path <- paste0(c(0, buffers), collapse = "_") %>%
    paste0(file.path(cache_folder, "gbif_thinned_cache_"), ., ".csv")
  thinned_cache <- get_gbif_thinned_cache(thinned_cache_path, stack_path, buffers)
  
  # Entities to thin:
  entities_to_thin <- thinned_cache[, .(genusKey, speciesKey, taxonKey)] %>%
    data.table::fsetdiff(entities, .) %>%
    unique()
  
  # Thin entities that need to be thinned:
  if (nrow(entities_to_thin) > 0) {
    # Create log file:
    log_path <- file.path(
      occ_folder, 
      paste0("gbif_thinning_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
    )
    if (file.exists(log_path)) file.remove(log_path)
    file.create(log_path, showWarnings = FALSE)
    log_file <- file(log_path, open = "wt")
    sink(file = log_file, type = "message")
    message(paste("LOG FILE - THINNING", nrow(entities_to_thin), 
                  "GBIF entities -", Sys.time()))
    sink(type =  "message")
    close(log_file)
    
    # Create a raster brick from the raster stack (to ensure faster retrieval):
    if (use_raster_brick) {
      message("Creating raster brick object from raster stack...")
      brick <- raster::brick(stack)
      message("Raster brick object successfully created!")
    } else {
      message(paste("Informational message: if enough memory is available,",
                    "consider setting use_raster_brick = TRUE.",
                    "Bioclimatic variables retrieval is faster from raster brick!"))
      brick <- stack
    }
    
    # Thin and retrieve bioclimatic variables:
    nb_cores <- parallel::detectCores()
    message(paste("Thinning and retrieving bioclimatic variables for", 
                  nrow(entities_to_thin), "GBIF entities..."))
    pb <- progress::progress_bar$new(
      total = nrow(entities_to_thin),
      format = "  :current/:total :percent |:bar| :elapsedfull",
      incomplete = " ",
      force = TRUE
    )
    parallel::mclapply(1:nrow(entities_to_thin), function(x) {
      tryCatch({
        thin_retrieve_gbif_single_entity(
          entities_to_thin[x, ], occ_path, occ_table, thinned_path,
          thinned_cache_path, brick, buffers
        )
      }, error = function(e) {
        cat("Entity genusKey=", entities_to_thin[x,][["genusKey"]],
            ", speciesKey=", entities_to_thin[x,][["speciesKey"]],
            ", taxonKey=", entities_to_thin[x,][["taxonKey"]],
            " issued an error!\n\t", conditionMessage(e), "\n",
            file = log_path, append = TRUE, sep = "")
      }, finally = {})
      if(x %% nb_cores == 0) pb$update(x/nrow(entities_to_thin))
    }, mc.cores = nb_cores)
    pb$terminate()
    message()
  }
  
  # Return last modification time of thinned occurrences database:
  file.info(thinned_path)$mtime
}

#' Get connection of the GBIF thinned entities database
#' 
#' This function returns a database connection to the GBIF thinned entities 
#' database. It checks that the database is up-to-date with stack files and 
#' (re)creates it if necessary.
#' 
get_gbif_thinned_db_con <- function(db_path, stack_path, brick) {
  # Create cache database if it does not exist or is outdated:
  if (file.exists(db_path)) {
    if (file.info(db_path)$ctime > file.info(stack_path)$ctime) {
      db <- RSQLite::dbConnect(RSQLite::SQLite(), db_path)
    } else {
      file.remove(db_path)
      db <- create_gbif_thinned_db(db_path, brick)
    }
  } else {
    db <- create_gbif_thinned_db(db_path, brick)
  }
  db
}

#' Create a database for thinned GBIF occurrences with bioclimatic variables
#' 
create_gbif_thinned_db <- function(
  db_path, brick, table_name = "thinned"
) {
  # Create cache folder if it does not exist:
  dir.create(dirname(db_path), showWarnings = FALSE, recursive = TRUE)
  
  # Create the database:
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = db_path)
  paste0(names(brick), collapse = ", ") %>%
    paste0("CREATE TABLE ", table_name, 
           " (genusKey INTEGER, speciesKey INTEGER, taxonKey INTEGER, ",
           "cell INTEGER, ", ., ")") %>%
    RSQLite::dbExecute(db, .)
  db
}

#' Return the cache of thinned GBIF entities or (re)create it if outdated
#' 
get_gbif_thinned_cache <- function(cache_path, stack_path, buffers) {
  if (file.exists(cache_path)) {
    if (file.info(cache_path)$ctime > file.info(stack_path)$ctime) {
      cache <- data.table::fread(cache_path, na.strings = c("", "NA"), 
                                 colClasses = "integer")
    } else {
      file.remove(cache_path)
      cache <- create_gbif_thinned_cache(cache_path, buffers)
    }
  } else {
    cache <- create_gbif_thinned_cache(cache_path, buffers)
  }
  cache
}

#' Create the cache of GBIF entities that have already been thinned
#' 
create_gbif_thinned_cache <- function(cache_path, buffers) {
  cache <- data.table::data.table(
    genusKey = integer(),
    speciesKey = integer(),
    taxonKey = integer(),
    nb_occurrences = integer(),
    nb_thinned = integer()
  )
  cache[, paste0("nb_buffer_", c(0, buffers)) := integer()]
  data.table::fwrite(cache, cache_path)
  cache
}


# Thin and retrieve bioclimatic conditions of a single GBIF entity =============
#' Thin and retrieve bioclimatic conditions of a single GBIF entity
#' 
thin_retrieve_gbif_single_entity <- function(
  entity, occ_path, occ_table, thinned_path, cache_path, brick, buffers
) {
  # Check if the entities database up-to-date and create it if necessary:
  thinned_db <- RSQLite::dbConnect(RSQLite::SQLite(), thinned_path)
  on.exit(RSQLite::dbDisconnect(thinned_db))
  
  # Open connection with GBIF occurrences database:
  occ_db <- RSQLite::dbConnect(RSQLite::SQLite(), occ_path)
  on.exit(RSQLite::dbDisconnect(occ_db), add = TRUE)
  
  # Get occurrences for the given entity:
  if (!is.na(entity[['speciesKey']])) {
    occurrences <- paste0(
      "SELECT DISTINCT decimalLongitude, decimalLatitude FROM ", occ_table,
      " WHERE genusKey = ", entity[['genusKey']], 
      " AND speciesKey = ", entity[['speciesKey']], 
      " AND taxonKey = ", entity[['taxonKey']]
    ) %>% 
      RSQLite::dbGetQuery(occ_db, .) %>% 
      data.table::as.data.table()
  } else {
    occurrences <- paste0(
      "SELECT DISTINCT decimalLongitude, decimalLatitude FROM ", occ_table,
      " WHERE genusKey = ", entity[['genusKey']], 
      " AND taxonKey = ", entity[['taxonKey']],
      " AND speciesKey IS NULL"
    ) %>% 
      RSQLite::dbGetQuery(occ_db, .) %>% 
      data.table::as.data.table()
  }
  nb_occurrences <- nrow(occurrences)
  
  # Proceed if the occurrences table is not empty:
  if (nb_occurrences > 0) {
    # Thin occurrences and retrieve their bioclimatic variables:
    occurrences %<>% get_thinned_bioclim_wo_cache(brick, buffers)
    nb_thinned <- nrow(occurrences)
    nb_buffers <- c(0, buffers) %>% 
      lapply(function (x) {sum(occurrences[['buffer']] == x)})
    
    # Append outcome to the thinned occurrences data.table:
    occurrences[, ':='(
      genusKey = entity[['genusKey']], 
      speciesKey = entity[['speciesKey']], 
      taxonKey = entity[['taxonKey']]
    )]
    RSQLite::dbAppendTable(thinned_db, "thinned", occurrences[, -"buffer"], 
                           "PRAGMA busy_timeout = 20 * 1000")
  } else {
    nb_thinned <- 0
    nb_buffers <- rep(0, times = length(buffers) + 1)
  }
  
  # Write statistics to cache:
  cache <- data.table::data.table(
    genusKey = entity[['genusKey']],
    speciesKey = entity[['speciesKey']],
    taxonKey = entity[['taxonKey']],
    nb_occurrences = nb_occurrences,
    nb_thinned = nb_thinned
  )
  cache[, paste0("nb_buffer_", c(0, buffers)) := nb_buffers]
  data.table::fwrite(cache, cache_path, append = TRUE)
}


# Retrieve Web of Life bioclimatic variables ===================================
retrieve_wol_bioclim <- function(
  wol_metadata, prob_nets, stack, use_raster_brick = TRUE, buffers = c(5000, 10000)
) {
  # Remove problematic networks:
  wol_metadata %<>% .[!is.na(lon) & !is.na(lat) & !net_name %in% prob_nets] %>%
    unique(by = "loc_id")
  
  # Order buffers:
  buffers %<>% unique() %>% sort()
  
  # Create a raster brick from the raster stack (to ensure faster retrieval):
  if (use_raster_brick) {
    message("Creating raster brick object from raster stack...")
    brick <- raster::brick(stack)
    message("Raster brick object successfully created!")
  } else {
    message(paste("Informational message: if enough memory is available,",
                  "consider setting use_raster_brick = TRUE.",
                  "Bioclimatic variables retrieval is faster from raster brick!"))
    brick <- stack
  }
  
  # Get bioclimatic conditions:
  wol_bioclim <- wol_metadata[, .(lon, lat)] %>%
    get_thinned_bioclim_wo_cache(brick, buffers)
  
  # Merge bioclimatic conditions with Web of Life locations:
  wol_metadata_w_grid <- wol_metadata[, .(
    loc_id, 
    cell = raster::cellFromXY(brick, data.table::data.table(lon, lat))
  )]
  merge(wol_metadata_w_grid, wol_bioclim[, -"buffer"], by = c("cell"), all.x = TRUE)
}
