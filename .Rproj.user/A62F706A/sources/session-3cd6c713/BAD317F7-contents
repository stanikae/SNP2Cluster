# First install tidyverse

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



# Source conf file --------------------------------------------------------

# source("./conf/KPN_Temb_BabyGERMS_config.R")
source(paste0("./conf/",conf_file))


# Read data into data frames ----------------------------------------------

var_order <- c(Main_var,Var_01,Var_02)

datesDF <- read_csv(dates_path, col_names = T) %>%
  dplyr::select(1,any_of(var_order))

names(datesDF)[1] <- "sampleID"
# names(datesDF)[names(datesDF) %in% Var_02]

val1 <- which(names(datesDF) %in% Var_02)
names(datesDF)[val1] <- "TakenDate"
Var_02 <- "TakenDate"


# Run main program --------------------------------------------------------
# snpco=20
# daysco=45
# 
# # Optionally a vector of SNP threshols can be provided together with matching
# # day intervals in a second vector 
# 
# # snpco=c(11,20,20,25,11)
# # daysco=c(60,14,60,45,14)

comparisons <- data.frame(snp=snpco,days=daysco)

for(i in 1:nrow(comparisons)){
  snpco=comparisons[i,1]
  daysco=comparisons[i,2]

  # source("~/GitHub/SNP2Cluster/conf/BabyGERMS_kpn_temb.R")

  

  # rm(list=ls())
  src_path <- file.path(MainDirPath)
  setwd(src_path)
  path_dir_main<-file.path(getwd())
  source("./EDA_snp_analysis_v2.R")
}
