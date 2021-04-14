# Download raw Web of Life data ================================================
#' List networks from a given interaction types available in Web of Life
#' 
#' @description 
#' List all the networks that correspond to the specified interaction type. 
#' The download date specified is not taken into account (changing it ensures 
#' that the data is downloaded again). Possible interaction types are:
#' * 3: plant-ant
#' * 5: pollination
#' * 6: seed dispersal
#' * 8: host-parasite
#' * 7: food webs
#' * 10: plant-herbivore
#' * 11: anemone-fish
#' * All
#' 
list_wol_networks <- function(interaction_type = "All", download_date = NA) {
  # Resolve interaction type name:
  interaction_id <- interaction_type
  if (is.character(interaction_type)) {
    if (interaction_type != "All") {
      interaction_id <- switch(
        interaction_type,
        "plant-ant" = 3,
        "pollination" = 5,
        "seed-dispersal" = 6,
        "host-parasite" = 8,
        "food-webs" = 7,
        "plant-herbivore" = 10,
        "anemone-fish" = 11,
        "All" = "All"
      )
    }
  }
  
  # Get list of data and parse it with JSON:
  "http://www.web-of-life.es/networkslist.php?type=" %>%
    paste0(interaction_id, "&data=All") %>%
    readLines(warn = FALSE) %>%
    glue::glue_collapse() %>%
    rjson::fromJSON()
}
