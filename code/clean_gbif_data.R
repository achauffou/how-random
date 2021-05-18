# Extract raw GBIF occurrences =================================================
#' Extract raw GBIF archives to a single csv file
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

