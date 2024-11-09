
# Set work directory ------------------------------------------------------


if (!require("tidyverse", quietly = TRUE)){
  install.packages("tidyverse")
  library(tidyverse)
}

# NOTES -------------------------------------------------------------------
# rm(list=ls())
# Get path to conf file and set it under the "Source conf file" section below
thisPath <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  if (length(grep("^-f$", cmdArgs)) > 0) {
    # R console option
    normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
  } else if (length(grep("^--file=", cmdArgs)) > 0) {
    # Rscript/R console option
    scriptPath <- normalizePath(dirname(sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)])))[1]
  } else if (Sys.getenv("RSTUDIO") == "1") {
    # RStudio
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else if (is.null(attr(stub, "srcref")) == FALSE) {
    # 'source'd via R console
    dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
  } else {
    stop("Cannot find file path")
  }
}
# function source: https://gist.github.com/jasonsychau/ff6bc78a33bf3fd1c6bd4fa78bbf42e7

MainDirPath <- thisPath()
# Clear env ---------------------------------------------------------------

src_path <- file.path(MainDirPath)
setwd(src_path)
