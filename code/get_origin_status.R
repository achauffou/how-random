# Find out the country and region of each site =================================
#' Get the country or region of each Web of Life site
#' 
get_sites_regions_codes <- function(metadata, gadm_countries, wgsrpd_l3) {
  metadata[!is.na(lat)][
    , .(lon = unique(lon), lat = unique(lat)), by = .(loc_id)][, ':='(
    GADM_code = get_gadm_code(lon, lat, gadm_countries),
    WGSRPD_code = get_wgsrpd_code(lon, lat, wgsrpd_l3)
  )]
}

#' Get GADM country code of given coordinates
#' 
get_gadm_code <- function(lon, lat, gadm_countries) {
  data.frame(lon = lon, lat = lat) %>%
    sp::SpatialPoints(proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>%
    sp::over(gadm_countries) %$%
    GID_0 %>%
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
