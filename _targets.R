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
read_YAML_config_targets <- list(
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
  tar_target(wol_fun_groups_info_path, config$wol_fun_groups_info_path),
  tar_target(wol_interaction_type, config$wol_interaction_type),
  tar_target(wol_supp_data_names_path, config$wol_supp_data_names_path),
  tar_target(worldclim_download_url, config$worldclim_download_url)
)


# Download raw data ============================================================
# Web of life data -----
download_web_of_life_data_targets <- list(
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
download_itis_data_targets <- tar_target(
  itis_raw_archive,
  download_from_url(itis_download_url, "data/raw/itis_sqlite.zip",
                    download_date),
  format = "file"
)

# Terrestrial ecoregions -----
download_ecoregions_data_targets <- tar_target(
  ecoregions_raw_archive,
  download_from_url(ecoregions_download_url, 
                    "data/raw/terrestrial_ecoregions.zip", download_date),
  format = "file"
)

# Climate data -----
download_climate_data_targets <- list(
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
download_raw_data_targets <- list(
  download_web_of_life_data_targets,
  download_itis_data_targets,
  download_ecoregions_data_targets,
  download_climate_data_targets
)


# Read raw data ================================================================
# Web of Life -----
read_web_of_life_data_targets <- list(
  tar_target(
    wol_raw_data,
    read_raw_wol_data(wol_raw_archive)
  ),
  tar_target(wol_raw_networks, wol_raw_data$networks),
  tar_target(wol_raw_metadata, wol_raw_data$metadata)
)

# ITIS databases -----
read_itis_data_targets <- list(
  tar_target(
    itis_raw_data,
    read_raw_itis_data(itis_raw_archive)
  )
)

# Read manually created data -----
read_manual_data_targets <- list(
  tar_target(
    wol_supp_data_names_file,
    wol_supp_data_names_path,
    format = "file"
  ),
  tar_target(
    wol_supp_data_names,
    readLines(wol_supp_data_names_file)
  ),
  tar_target(
    wol_fun_groups_info_file,
    wol_fun_groups_info_path,
    format = "file"
  ),
  tar_target(
    wol_fun_groups_info,
    data.table::fread(wol_fun_groups_info_file)
  )
)

# List all read targets -----
read_raw_data_targets <- list(
  read_web_of_life_data_targets,
  read_itis_data_targets,
  read_manual_data_targets
)


# Prepare interactions data ====================================================
# Prepare metadata -----
prepare_interactions_metadata_targets <- list(
  tar_target(
    wol_metadata,
    create_wol_metadata_loc_id(wol_raw_metadata)
  )
)

# Clean species names -----
clean_species_names_targets <- list(
  tar_target(
    wol_raw_species,
    get_raw_wol_species(wol_networks_wo_supp_data)
  ),
  tar_target(
    wol_raw_species_w_metadata,
    add_metadata_to_wol_species(wol_raw_species, wol_metadata, wol_fun_groups_info)
  )
)

# Prepare interactions -----
prepare_interactions_targets <- list(
  tar_target(
    wol_networks_wo_supp_data,
    remove_supp_data_from_wol_networks(wol_raw_networks, wol_supp_data_names)
  )
)

# List all targets to prepare interactions data -----
prepare_interactions_data_targets <- list(
  prepare_interactions_metadata_targets,
  clean_species_names_targets,
  prepare_interactions_targets
)


# Compile TeX manuscripts ======================================================
compile_TeX_manuscripts_targets <-list(
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
  read_YAML_config_targets,
  download_raw_data_targets,
  read_raw_data_targets,
  prepare_interactions_data_targets,
  compile_TeX_manuscripts_targets
)
