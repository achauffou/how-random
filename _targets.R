# Configuration ================================================================
# Load targets package -----
library(targets)

# Load other crucial packages -----
library(magrittr)

# Set R options -----
QUIET_DOWNLOADS <- FALSE
options(download.file.method = "curl")
options(download.file.extra = "-L")
options(tinytex.engine = "lualatex")
options(tinytex.engine_args = "-shell-escape")
options(tinytex.bib_engine = "biber")
options(tinytex.compile.min_times = 2)

# Set options to load R packages listed in .r_packages if necessary -----
tar_option_set(packages = readLines(".r_packages"))

# Load all functions scripts in the code folder -----
f <- lapply(list.files("code", recursive = TRUE, full.names = TRUE), source)


# Read YAML configuration ======================================================
read_YAML_config <- list(
  tar_target(
    config_file,
    "config.yaml",
    format = "file"
  ),
  tar_target(
    config,
    yaml::read_yaml(config_file)
  ),
  tar_target(download_date, config$download_date),
  tar_target(ecoregions_download_url, config$ecoregions_download_url),
  tar_target(envirem_bioclim_download_url, config$envirem_bioclim_download_url),
  tar_target(envirem_topo_download_url, config$envirem_topo_download_url),
  tar_target(itis_download_url, config$itis_download_url),
  tar_target(tex_folders_to_compile, config$tex_folders_to_compile),
  tar_target(wol_interaction_type, config$wol_interaction_type),
  tar_target(worldclim_download_url, config$worldclim_download_url)
)


# Download raw data ============================================================
# Web of life data -----
download_web_of_life_data <- list(
  tar_target(
    wol_download_list,
    download_wol_networks_list(wol_interaction_type, download_date)
  ),
  tar_target(
    wol_raw_archive,
    download_wol_networks_raw_archive(
      wol_download_list,  "data/raw/wol.zip"
    ),
    format = "file"
  )
)

# ITIS database -----
download_itis_data <- tar_target(
  itis_raw_archive,
  download_from_url(itis_download_url, "data/raw/itis_sqlite.zip",
                    download_date),
  format = "file"
)

# Terrestrial ecoregions -----
download_ecoregions_data <- tar_target(
  ecoregions_raw_archive,
  download_from_url(ecoregions_download_url, 
                    "data/raw/terrestrial_ecoregions.zip", download_date),
  format = "file"
)

# Climate data -----
download_climate_data <- list(
  tar_target(
    worldclim_raw_archive,
    download_from_url(worldclim_download_url, 
                     "data/raw/worldclim_2-5.zip", download_date),
    format = "file"
  ),
  tar_target(
    envirem_bioclim_raw_archive,
    download_from_url(envirem_bioclim_download_url, 
                      "data/raw/envirem_bioclim_2-5.zip", download_date),
    format = "file"
  ),
  tar_target(
    envirem_topo_raw_archive,
    download_from_url(envirem_topo_download_url, 
                      "data/raw/envirem_topo_2-5.zip", download_date),
    format = "file"
  )
)

# List all download targets -----
download_raw_data <- list(
  download_web_of_life_data,
  download_itis_data,
  download_ecoregions_data,
  download_climate_data
)


# Compile TeX manuscripts ======================================================
compile_TeX_manuscripts <-list(
  tar_target(
    tex_source_files_to_watch,
    list.files(
      path = tex_folders_to_compile,
      pattern = c(".*[tex|bib]$"),
      recursive = TRUE,
      full.names = TRUE
    ),
    format = "file",
    pattern = map(tex_folders_to_compile)
  ),
  tar_target(
    compile_tex_manuscripts,
    latexmk_from_path(
      grep("main.tex$", tex_source_files_to_watch, value = TRUE)[1]
    ),
    pattern = map(tex_source_files_to_watch)
  )
)


# List all project targets to make =============================================
list(
  read_YAML_config,
  download_raw_data,
  compile_TeX_manuscripts
)
