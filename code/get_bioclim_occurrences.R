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
  if (file.exists(save_path)) {
    raster_files <- list.files(temp_dir, full.names = T)
    if (min(file.info(archives)$ctime) < max(file.info(raster_files)$mtime)) {
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
