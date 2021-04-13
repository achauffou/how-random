# Configuration ================================================================
# Load targets package -----
library(targets)

# Set environment variables and options -----
Sys.setenv(CMDSTAN="/usr/local/cmdstan")
options(tinytex.engine = "lualatex")
options(tinytex.engine_args = "-shell-escape")
options(tinytex.bib_engine = "biber")
options(tinytex.compile.min_times = 2)

# Load all functions scripts in the code folder -----
f <- lapply(list.files("code", recursive = TRUE, full.names = TRUE), source)


# Compile TeX manuscripts ======================================================
list(
  tar_target(
    tex_folders_to_compile,
    c("manuscript/outline/2021-04-07")
  ),
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
