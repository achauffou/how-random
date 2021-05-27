# Utilities to save objects to a file and return the object ====================
#' Save object to a file matching its class and return the object
#' 
save_obj <- function(obj, name, folder) {
  UseMethod("save_obj")
}

save_obj.default <- function(obj, name, folder) {
  # Create folder if it does not exist yet:
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  
  # Save the object as an RDS file:
  saveRDS(obj, file = paste0(folder, "/", name, ".rds"))
  
  # Return the data.table:
  obj
}

save_obj.data.table <- function(obj, name, folder) {
  # Create folder if it does not exist yet:
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  
  # Write the data.table to a csv:
  data.table::fwrite(obj, file = paste0(folder, "/", name, ".csv"))
  
  # Return the data.table:
  obj
}
