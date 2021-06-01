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
