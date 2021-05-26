# Configuration ================================================================
# Load targets package:
library(targets)

# Load other crucial packages:
library(magrittr)

# Set R options:
QUIET_DOWNLOADS <- FALSE
options(download.file.method = "curl")
options(download.file.extra = "-L")
options(tinytex.engine = "lualatex")
options(tinytex.engine_args = "-shell-escape")
options(tinytex.bib_engine = "biber")
options(tinytex.compile.min_times = 2)
options(CHUNK_SIZE = 2E4)

# Set options to load R packages listed in .r_packages if necessary:
tar_option_set(packages = readLines(".r_packages"))

# Load all functions scripts in the code folder:
f <- list.files("code", recursive = TRUE, full.names = TRUE) %>%
  grep("\\.R$", ., value = TRUE) %>%
  lapply(source)


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
  tar_target(accepted_ranks, config$accepted_ranks),
  tar_target(aggregation_level, config$aggregation_level),
  tar_target(download_date, config$download_date),
  tar_target(ecoregions_download_url, config$ecoregions_download_url),
  tar_target(envirem_bioclim_download_url, config$envirem_bioclim_download_url),
  tar_target(envirem_topo_download_url, config$envirem_topo_download_url),
  tar_target(fun_groups_plausible_kingdoms, config$fun_groups_plausible_kingdoms),
  tar_target(itis_download_url, config$itis_download_url),
  tar_target(min_locations_per_species, config$min_locations_per_species),
  tar_target(stan_sim_specs, config$stan_simulations_specs),
  tar_target(tex_folders_to_compile, config$tex_folders_to_compile),
  tar_target(wol_aquatic_networks, config$wol_aquatic_networks),
  tar_target(wol_interaction_type, config$wol_interaction_type),
  tar_target(worldclim_download_url, config$worldclim_download_url)
)


# Download raw data ============================================================
# Web of life data:
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

# ITIS database:
download_itis_data_targets <- tar_target(
  itis_raw_archive,
  download_from_url(itis_download_url, "data/raw/itis_sqlite.zip",
                    download_date),
  format = "file"
)

# Terrestrial ecoregions:
download_ecoregions_data_targets <- tar_target(
  ecoregions_raw_archive,
  download_from_url(ecoregions_download_url,
                    "data/raw/terrestrial_ecoregions.zip", download_date),
  format = "file"
)

# Climate data:
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

# Rnaturalearth data:
download_rnaturalearth_targets <- list(
  tar_target(
    rnaturalearth_land_data_download,
    rnaturalearth::ne_download(
      type = "land", category = "physical", returnclass = "sp", scale = 10,
      load = FALSE, destdir = "data/raw/rnaturalearth"
    ) %>% list(file_name = ., download_date = download_date)
  )
)

# GBIF occurrence data:
download_occurrence_data_targets <- list(
  tar_target(
    gbif_names_to_suggest,
    select_names_to_gbif_suggest(wol_species, taxonomic_dict,
                                 min_locations_per_species, aggregation_level)
  ),
  tar_target(
    gbif_keys_dict,
    suggest_gbif_names(gbif_names_to_suggest, "data/cache/gbif_keys_cache.csv")
  ),
  tar_target(
    gbif_keys,
    select_gbif_keys_to_download(gbif_names_to_suggest, gbif_keys_dict,
                                 accepted_ranks)
  ),
  tar_target(
    gbif_raw_archives,
    get_gbif_occurrences(gbif_keys[['gbif_key']], "data/raw/gbif",
                         "data/cache/gbif_downloads_cache.csv", download_date)
  )
)

# List all download targets:
download_raw_data_targets <- list(
  download_web_of_life_data_targets,
  download_itis_data_targets,
  download_ecoregions_data_targets,
  download_climate_data_targets,
  download_rnaturalearth_targets,
  download_occurrence_data_targets
)


# Read raw data ================================================================
# Web of Life:
read_web_of_life_data_targets <- list(
  tar_target(
    wol_raw_data,
    read_raw_wol_data(wol_raw_archive)
  ),
  tar_target(wol_raw_networks, wol_raw_data$networks),
  tar_target(wol_raw_metadata, wol_raw_data$metadata)
)

# ITIS databases:
read_itis_data_targets <- list(
  tar_target(
    itis_raw_data,
    read_raw_itis_data(itis_raw_archive)
  )
)

# Rnaturalearth data:
read_rnaturalearth_targets <- list(
  tar_target(
    rnaturalearth_land_data,
    rnaturalearth::ne_load(
      destdir = "data/raw/rnaturalearth",
      file_name = rnaturalearth_land_data_download[['file_name']]
    )
  )
)

# GBIF occurrences:
read_gbif_data_targets <- list(
  tar_target(
    gbif_last_cleaning_update,
    process_gbif_raw_archives(
      gbif_raw_archives,
      "data/processed/gbif",
      rnaturalearth_land_data,
      "data/cache/gbif_processing_cache.csv"
    )
  )
)

# Read manually created data:
read_manual_data_targets <- list(
  tar_target(
    wol_supp_data_names_file,
    "data/wol_supp_data_names.txt",
    format = "file"
  ),
  tar_target(
    wol_supp_data_names,
    readLines(wol_supp_data_names_file)
  ),
  tar_target(
    wol_fun_groups_info_file,
    "data/wol_fun_groups_info.csv",
    format = "file"
  ),
  tar_target(
    wol_fun_groups_info,
    data.table::fread(wol_fun_groups_info_file, na.strings = c("", "NA"))
  ),
  tar_target(
    wol_manual_locations_file,
    "data/wol_manual_locations.csv",
    format = "file"
  ),
  tar_target(
    wol_manual_locations,
    data.table::fread(wol_manual_locations_file, na.strings = c("", "NA"))
  ),
  tar_target(
    wol_manual_species_names_file,
    "data/wol_manual_species_names.csv",
    format = "file"
  ),
  tar_target(
    wol_manual_species_names,
    data.table::fread(wol_manual_species_names_file, na.strings = c("", "NA"))
  )
)

# List all read targets:
read_raw_data_targets <- list(
  read_web_of_life_data_targets,
  read_itis_data_targets,
  read_rnaturalearth_targets,
  read_gbif_data_targets,
  read_manual_data_targets
)


# Prepare interactions data ====================================================
# Prepare metadata:
prepare_interactions_metadata_targets <- list(
  tar_target(
    wol_metadata,
    create_wol_metadata_loc_id(wol_raw_metadata, wol_manual_locations)
  )
)

# Clean species names:
clean_species_names_targets <- list(
  tar_target(
    wol_raw_species,
    get_raw_wol_species(wol_networks_wo_supp_data)
  ),
  tar_target(
    wol_raw_species_w_metadata,
    add_metadata_to_wol_species(wol_raw_species, wol_metadata, wol_fun_groups_info)
  ),
  tar_target(
    wol_proposed_names,
    prepare_names_to_verify(wol_raw_species, wol_manual_species_names)
  ),
  tar_target(
    taxonomic_dict,
    check_proposed_names(wol_proposed_names, itis_raw_data,
                         "data/cache/taxonomic_dict_cache.csv")
  ),
  tar_target(
    wol_verified_names,
    get_verified_names(wol_proposed_names, taxonomic_dict)
  ),
  tar_target(
    wol_species_cleaned,
    select_verified_species(
      wol_raw_species_w_metadata, wol_verified_names,
      fun_groups_plausible_kingdoms, accepted_ranks, aggregation_level
    )
  ),
  tar_target(
    wol_problematic_networks,
    detect_problematic_networks(wol_species_cleaned, wol_metadata) %>%
      c(wol_aquatic_networks)
  ),
  tar_target(
    wol_species,
    remove_problematic_species(wol_species_cleaned, wol_problematic_networks)
  )
)

# Prepare interactions:
prepare_interactions_targets <- list(
  tar_target(
    wol_networks_wo_supp_data,
    remove_supp_data_from_wol_networks(wol_raw_networks, wol_supp_data_names)
  ),
  tar_target(
    wol_interactions,
    get_wol_interactions(wol_networks_wo_supp_data, wol_metadata,
                         wol_species, wol_fun_groups_info,
                         wol_problematic_networks)
  )
)

# List all targets to prepare interactions data:
prepare_interactions_data_targets <- list(
  prepare_interactions_metadata_targets,
  clean_species_names_targets,
  prepare_interactions_targets
)


# Simulate Stan models =========================================================
simulate_stan_models_targets <- list(
  tar_target(
    stan_sim_data,
    generate_stan_sim_data(stan_sim_specs[[1]], "results/simulations"),
    pattern = map(stan_sim_specs),
    format = "file"
  ),
  tar_target(
    stan_sim_starts,
    generate_stan_sim_start_values(stan_sim_specs[[1]], "results/simulations"),
    pattern = map(stan_sim_specs),
    format = "file"
  )
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
  simulate_stan_models_targets,
  compile_TeX_manuscripts_targets
)
