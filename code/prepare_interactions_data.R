# Prepare interactions metadata ================================================
#' Add a column with location ID to Web of Life metadata
#' 
create_wol_metadata_loc_id <- function(metadata) {
  metadata[order(net_name), ':='(loc_id = paste0("wol_", .GRP)), by = .(lat, lon)]
}
