# Read Web of Life data ========================================================
#' Read a raw Web of Life archive and return its networks and metadata
#' 
read_raw_wol_data <- function(file_path, format = "csv") {
  # Unzip the raw archive:
  dir_path <- paste0(tempdir(), "/", 
                     tools::file_path_sans_ext(basename(file_path)))
  unzip(file_path, exdir = dir_path)
  
  # Return invisibly a list with networks and metadata:
  list(
    networks = read_raw_wol_networks(dir_path, format = "csv"),
    metadata = read_raw_wol_metadata(
      file.path(dir_path, paste0("references.", format))
    )
  )
}

#' Read raw Web of Life networks located in an unzipped directory
#' 
read_raw_wol_networks <- function(dir_path, format = "csv") {
  # Get file and network names of all networks in the folder:
  regexp_pattern <- paste0("^(?!(references|README))(\\w*)(\\.", format, ")$")
  file_names <- list.files(dir_path) %>%
    grep(regexp_pattern, ., perl = TRUE, value = TRUE)
  networks_names <- tools::file_path_sans_ext(file_names)
  
  # Read network files:
  suppressMessages({
    suppressWarnings({
      file_names %>% 
        file.path(dir_path, .) %>%
        purrr::map(readr::read_csv) %>%
        purrr::map(~dplyr::filter(., X1 != 'Abundance"')) %>%
        purrr::map(format_wol_network_as_matrix) %>%
        `names<-`(networks_names)
    })
  })
}

#' Turn an interaction network from data.frame to matrix
#' 
format_wol_network_as_matrix <- function(network) {
  column_names <- names(network)[-1]
  row_names <- as.character(network$X1)
  as.matrix(network[, -1]) %>%
    `rownames<-`(row_names) %>%
    `colnames<-`(column_names)
}

#' Read raw Web of Life metadata located in a csv references file
#' 
read_raw_wol_metadata <- function(file_path) {
  suppressMessages({
    readr::read_csv(
      file_path,
      col_names = c("net_name", "n_spp", "n_int", "c", "int_type", "data_type", 
        "reference", "loc_name", "lat", "lon"),
      skip = 1
    ) %>%
    data.table::as.data.table()
  })
}


# Read ITIS database ===========================================================
#' Read raw ITIS database archive
#' 
#' Read an ITIS sqlite database archive and return taxonomic units and synonym 
#' links datatables.
#' 
read_raw_itis_data <- function(file_path) {
  # Unzip the raw archive:
  dir_path <- paste0(tempdir(), "/", 
                     tools::file_path_sans_ext(basename(file_path)))
  unzipped_files <- unzip(file_path, exdir = dir_path)
  
  # Establish connection with the database:
  db_path <- unzipped_files[stringr::str_detect(unzipped_files, "ITIS")]
  db_con <- RSQLite::dbConnect(RSQLite::SQLite(), db_path)
  
  # Retrieve table with taxonomic units:
  taxonomic_units <- RSQLite::dbReadTable(db_con, "taxonomic_units") %>%
    data.table::as.data.table() %>%
    # Keep only plants, animals, and drop useless columns:
    .[kingdom_id %in% c(3,5), .(tsn, complete_name, n_usage, rank_id, 
                                unaccept_reason)]
  # Retrieve table with synonym links:
  synonym_links <- RSQLite::dbReadTable(db_con, "synonym_links") %>%
    data.table::as.data.table()
  
  # Close connection and return tables:
  RSQLite::dbDisconnect(db_con)
  list(taxonomic_units = taxonomic_units, synonym_links = synonym_links)
}
