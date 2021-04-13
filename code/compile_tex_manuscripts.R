#' Call tinytex::latexmk() from within the file directory
#' 
#' Call tinytex::latexmk() from within the directory of a TeX manuscript main 
#' source file to compile the it.
#' 
latexmk_from_path <- function(file_path) {
  wd <- getwd()
  setwd(dirname(file_path))
  out <- tinytex::latexmk(basename(file_path))
  setwd(wd)
  out
}
