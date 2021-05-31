# Extract and stack bioclimatic rasters ========================================
#' Extract and stack together raster files from multiple archives
#' 
stack_bioclim_archives <- function(
  archives, extent_coords, save_path, temp_dir = file.path(tempdir(), "bioclim")
) {
  # Create output folder if it does not exist:
  dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
  
  # Skip the task if the rasters have already been extracted and are up-to-date:
  if (file.exists(save_path)) {
    if (min(file.info(archives)$ctime) < file.info(save_path)$ctime) {
      return(save_path)
    }
  }
  
  # Extract archives to temporary directory:
  lapply(archives, unzip, exdir = temp_dir)
  
  # List all raster files:
  raster_files <- list.files(temp_dir, full.names = T)
  
  # Stack together all rasters and correct their extents:
  list.files(temp_dir, full.names = T) %>% 
    lapply(crop_raster_from_file, extent_coords) %>%
    raster::stack() %>% 
    saveRDS(save_path)
  save_path
}

#' Read a raster file and crop it to a given extent
#' 
crop_raster_from_file <- function(file_path, extent_coords) {
  out <- raster::stack(file_path)
  ext <- as(raster::extent(extent_coords), 'SpatialPolygons')
  raster::crs(ext) <- raster::crs(out)
  raster::crop(out, ext)
}
