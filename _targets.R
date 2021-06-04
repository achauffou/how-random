# Configuration ================================================================
# Load targets package:
library(targets)

# Load other crucial packages:
library(magrittr)

# Set R options:
if (is.null(getOption("download.file.method"))) {
  options(download.file.method = "curl")
  options(download.file.extra = "-L")
}
if (is.null(getOption("tinytex.engine"))) {
  options(tinytex.engine = "lualatex")
  options(tinytex.engine_args = "-shell-escape")
}
options(tinytex.bib_engine = getOption("tinytex.bib_engine", default = "biber"))
options(tinytex.compile.min_times = getOption("tinytex.compile.min_times", default = 2))

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
  tar_target(bioclim_extent, config$bioclim_extent),
  tar_target(download_date, config$download_date),
  tar_target(ecoregions_download_url, config$ecoregions_download_url),
  tar_target(envirem_bioclim_download_url, config$envirem_bioclim_download_url),
  tar_target(envirem_topo_download_url, config$envirem_topo_download_url),
  tar_target(manual_data_folder, config$folder_structure$manual_data),
  tar_target(cache_folder, config$folder_structure$cache),
  tar_target(raw_data_folder, config$folder_structure$raw_data),
  tar_target(processed_data_folder, config$folder_structure$processed_data),
  tar_target(stan_src_folder, config$folder_structure$stan_sources),
  tar_target(stan_bin_folder, config$folder_structure$stan_binaries),
  tar_target(results_sim_folder, config$folder_structure$results_simulations),
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
      wol_download_list, file.path(raw_data_folder, "wol.zip"), 
      download_date = download_date
    ),
    format = "file"
  )
)

# ITIS database:
download_itis_data_targets <- tar_target(
  itis_raw_archive,
  download_from_url(
    itis_download_url, file.path(raw_data_folder, "itis_sqlite.zip"), 
    download_date
  ),
  format = "file"
)

# Terrestrial ecoregions:
download_ecoregions_data_targets <- tar_target(
  ecoregions_raw_archive,
  download_from_url(
    ecoregions_download_url, 
    file.path(raw_data_folder, "terrestrial_ecoregions.zip"), download_date),
  format = "file"
)

# Climate data:
download_climate_data_targets <- list(
  tar_target(
    worldclim_raw_archive,
    download_from_url(
      worldclim_download_url, file.path(raw_data_folder, "worldclim_2-5.zip"), 
      download_date
    ),
    format = "file"
  ),
  tar_target(
    envirem_bioclim_raw_archive,
    download_from_url(
      envirem_bioclim_download_url, 
      file.path(raw_data_folder, "envirem_bioclim_2-5.zip"), download_date
    ),
    format = "file"
  ),
  tar_target(
    envirem_topo_raw_archive,
    download_from_url(
      envirem_topo_download_url,
      file.path(raw_data_folder, "envirem_topo_2-5.zip"), download_date),
    format = "file"
  )
)

# Rnaturalearth data:
download_rnaturalearth_targets <- list(
  tar_target(
    rnaturalearth_land_data,
    download_rnaturalearth_land_data(raw_data_folder, 10, download_date),
    format = "file"
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
    suggest_gbif_names(
      gbif_names_to_suggest, 
      file.path(cache_folder, "gbif_keys_cache.csv")
    )
  ),
  tar_target(
    gbif_keys,
    select_gbif_keys_to_download(gbif_names_to_suggest, gbif_keys_dict,
                                 accepted_ranks) %>%
      save_obj("gbif_keys", processed_data_folder)
  ),
  tar_target(
    gbif_raw_archives,
    get_gbif_occurrences(
      gbif_keys[['gbif_key']], 
      file.path(raw_data_folder, "gbif"),
      file.path(cache_folder, "gbif_downloads_cache.csv"), 
      download_date
    )
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

# GBIF occurrences:
read_gbif_data_targets <- list(
  tar_target(
    gbif_last_cleaning_update,
    process_gbif_raw_archives(
      gbif_raw_archives,
      file.path(processed_data_folder, "gbif"),
      readRDS(rnaturalearth_land_data),
      file.path(cache_folder, "gbif_processing_cache.csv")
    )
  )
)

# Read manually created data:
read_manual_data_targets <- list(
  tar_target(
    wol_supp_data_names_file,
    file.path(manual_data_folder, "wol_supp_data_names.txt"),
    format = "file"
  ),
  tar_target(
    wol_supp_data_names,
    readLines(wol_supp_data_names_file)
  ),
  tar_target(
    wol_fun_groups_info_file,
    file.path(manual_data_folder, "wol_fun_groups_info.csv"),
    format = "file"
  ),
  tar_target(
    wol_fun_groups_info,
    data.table::fread(wol_fun_groups_info_file, na.strings = c("", "NA"))
  ),
  tar_target(
    wol_manual_locations_file,
    file.path(manual_data_folder, "wol_manual_locations.csv"),
    format = "file"
  ),
  tar_target(
    wol_manual_locations,
    data.table::fread(wol_manual_locations_file, na.strings = c("", "NA"))
  ),
  tar_target(
    wol_manual_species_names_file,
    file.path(manual_data_folder, "wol_manual_species_names.csv"),
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
                         file.path(cache_folder, "taxonomic_dict_cache.csv"))
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
    wol_species_cleaned %>%
      remove_problematic_species(wol_problematic_networks) %>%
      get_species_avg_rel_degree() %>%
      save_obj("wol_species", processed_data_folder)
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
                         wol_problematic_networks) %>%
      save_obj("wol_interactions", processed_data_folder)
  )
)

# List all targets to prepare interactions data:
prepare_interactions_data_targets <- list(
  prepare_interactions_metadata_targets,
  clean_species_names_targets,
  prepare_interactions_targets
)


# Perform bioclimatic suitability analyses =====================================
# Stack climatic data:
get_bioclim_stack_targets <- list(
  tar_target(
    raw_bioclim_stack,
    stack_bioclim_archives(
      c(worldclim_raw_archive, envirem_bioclim_raw_archive, 
        envirem_topo_raw_archive), 
      bioclim_extent, 
      file.path(processed_data_folder, "bioclim_vars"),
      download_date = download_date
    )
  )
)

# Thin and retrieve bioclimatic conditions:
thin_retrieve_gbif_bioclim_targets <- list(
  tar_target(
    gbif_entities_to_thin,
    get_gbif_entities_to_thin(
      gbif_keys, file.path(processed_data_folder, "gbif"), 
      last_occ_update = gbif_last_cleaning_update
    )
  ),
  tar_target(
    last_gbif_bioclim_update,
    thin_retrieve_gbif_entities(
      gbif_entities_to_thin, file.path(processed_data_folder, "gbif"),
      raw_bioclim_stack, cache_folder,
      file.path(processed_data_folder, "bioclim_vars")
    )
  ),
  tar_target(
    wol_bioclim,
    retrieve_wol_bioclim(wol_metadata, wol_problematic_networks, raw_bioclim_stack)
  ),
  tar_target(
    nb_occurrences_per_species,
    count_nb_occs_per_species(
      gbif_keys, wol_species, wol_bioclim, 
      file.path(processed_data_folder, "gbif", "bioclim.sqlite"), 
      file.path(cache_folder, "nb_occs_per_species_cache.csv"),
      aggregation_level = aggregation_level
    )
  )
)

# List all targets to perform bioclimatic suitability analyses:
perform_bioclim_analyses_targets <- list(
  get_bioclim_stack_targets,
  thin_retrieve_gbif_bioclim_targets
)


# Simulate Stan models =========================================================
simulate_stan_models_targets <- list(
  tar_target(
    stan_sim_data,
    generate_stan_sim_data(stan_sim_specs[[1]], results_sim_folder),
    pattern = map(stan_sim_specs),
    format = "file"
  ),
  tar_target(
    stan_sim_starts,
    generate_stan_sim_start_values(stan_sim_specs[[1]], results_sim_folder),
    pattern = map(stan_sim_specs),
    format = "file"
  ),
  tar_target(
    stan_sim_sources,
    stan_sim_specs[[1]]$stan_model %>%
      paste0(".stan") %>%
      file.path(stan_src_folder, .),
    pattern = map(stan_sim_specs),
    format = "file"
  ),
  tar_target(
    stan_sim_fits,
    run_stan_model(
      stan_sim_specs[[1]],
      readRDS(stan_sim_data),
      readRDS(stan_sim_starts),
      stan_sim_sources,
      stan_bin_folder,
      results_sim_folder
    ),
    pattern = map(stan_sim_specs, stan_sim_data, stan_sim_starts, stan_sim_sources),
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
  perform_bioclim_analyses_targets,
  simulate_stan_models_targets,
  compile_TeX_manuscripts_targets
)
