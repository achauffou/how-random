# Configuration ================================================================
# Load targets package -----
library(targets)

# Load other crucial packages -----
library(magrittr)

# Set environment variables and options -----
Sys.setenv(CMDSTAN="/usr/local/cmdstan")
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
  tar_target(tex_folders_to_compile, config$tex_folders_to_compile),
  tar_target(download_date, config$download_date),
  tar_target(wol_interaction_type, config$wol_interaction_type)
)


# Download raw data ============================================================
# Web of life data -----
download_web_of_life_data <- list(
  tar_target(
    wol_networks_list, 
    list_wol_networks(wol_interaction_type, download_date)
  )
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
    compile_outline,
    latexmk_from_path(
      grep("main.tex$", tex_source_files_to_watch, value = TRUE)[1]
    ),
    pattern = map(tex_source_files_to_watch)
  )
)


# List all targets to make =====================================================
list(
  read_YAML_config, 
  download_web_of_life_data,
  compile_TeX_manuscripts
)
