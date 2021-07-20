# Find out the country and region of each site =================================
#' Get the country or region of each Web of Life site
#' 
get_sites_regions_codes <- function(metadata, wab_countries, wgsrpd_l3) {
  metadata[!is.na(lat)][
    , .(lon = unique(lon), lat = unique(lat)), by = .(loc_id)][, ':='(
      WAB_code = get_wab_code(lon, lat, wab_countries),
      WGSRPD_code = get_wgsrpd_code(lon, lat, wgsrpd_l3)
    )]
}

#' Get WAB country code of given coordinates
#' 
get_wab_code <- function(lon, lat, wab_countries) {
  data.frame(lon = lon, lat = lat) %>%
    sp::SpatialPoints(proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
    sp::over(wab_countries) %$%
    iso3 %>%
    countrycode::countrycode(origin = "iso3c", destination = "iso2c")
}

#' Get WGSRPD region code of given coordinates
#' 
get_wgsrpd_code <- function(lon, lat, wgsrpd_l3) {
  data.frame(lon = lon, lat = lat) %>%
    sp::SpatialPoints(proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
    sp::over(wgsrpd_l3) %$%
    LEVEL3_COD
}
