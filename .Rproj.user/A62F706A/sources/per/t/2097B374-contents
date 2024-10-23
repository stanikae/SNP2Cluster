library(tidyverse)

# NOTES -------------------------------------------------------------------
rm(list=ls())
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

# source("~/GitHub/SNP2Cluster/conf/BabyGERMS_kpn_temb.R")
# source("~/GitHub/SNP2Cluster/conf/Thabo_p_aerogenosa.R")
# source("~/GitHub/SNP2Cluster/conf/Richael_s_pneumo.R")
source("./conf/KPN_Temb_BabyGERMS_config.R")


# Read data into data frames ----------------------------------------------
# # Thabo
# var_order <- c("sampleID",Main_var,Var_01)
# 
# var_sel <- c("sampleid","biosample","collection_year","collection_date","Country")
# mx <- read_excel(dates_path) %>% dplyr::filter(Host=="Human") %>% ncol()
# datesDF <- read_excel(dates_path) %>% 
#   dplyr::filter(Host=="Human") %>% 
#   dplyr::select(any_of(var_sel)) #%>%
#   # dplyr::select(any_of(var_order))
# # dplyr::filter(collection_date != "Missing")
# 
# names(datesDF)[1] <- "sampleID"


# get BabyGERMS metadata
# metdata now included in the dates file: 2023-05-01
mx <- read_csv(dates_path, col_names = F) %>% ncol()
datesDF <- read_csv(dates_path, col_names = F) %>% dplyr::select(all_of(mx),4,1,3,2) #,1,2,3,4)
names(datesDF) <- c("sampleID","TakenDate","Hospital","WardType","Ward")
datesDF$Ward <- dplyr::coalesce(datesDF$Ward,"UNKNOWN")
datesDF$WardType <- dplyr::coalesce(datesDF$WardType,"UNKNOWN")

names(datesDF)[1] <- "sampleID"




# # # Richael
# var_order <- c(Main_var,Var_01,Var_02)
# 
# datesDF <- read_csv(dates_path, col_names = T) %>%
#   dplyr::select(1,any_of(var_order))
# 
# names(datesDF)[1] <- "sampleID"
# # names(datesDF)[names(datesDF) %in% Var_02]
# 
# val1 <- which(names(datesDF) %in% Var_02)
# names(datesDF)[val1] <- "TakenDate"
# Var_02 <- "TakenDate"


# Run main program --------------------------------------------------------

# Richael-Ghana

# # Community-based analysis
# comparisons <- data.frame(snp=c(20),days=c(45))
# for(i in 1:nrow(comparisons)){
#   snpco=comparisons[i,1]
#   daysco=comparisons[i,2]
#   
#   # # Richael
#   source("~/GitHub/SNP2Cluster/conf/Richael_s_pneumo_community.R")
#   var_order <- c(Main_var,Var_01,Var_02)
#   
#   datesDF <- read_csv(dates_path, col_names = T) %>%
#     dplyr::select(1,any_of(var_order)) 
#   
#   names(datesDF)[1] <- "sampleID"
#   # names(datesDF)[names(datesDF) %in% Var_02]
#   
#   val1 <- which(names(datesDF) %in% Var_02)
#   names(datesDF)[val1] <- "TakenDate"
#   Var_02 <- "TakenDate"
#   
#   source("~/GitHub/SNP2Cluster/EDA_snp_analysis_v2.R")
#   
#   # rm(list=ls())
#   src_path <- file.path("C:/Users/Stanfordk/Documents/GitHub/SNP2Cluster")
#   setwd(src_path)
# }

# # Facility level analysis
# comparisons <- data.frame(snp=c(20),days=c(45))
# for(i in 1:nrow(comparisons)){
#   snpco=comparisons[i,1]
#   daysco=comparisons[i,2]
#   
#   # # Richael
#   source("~/GitHub/SNP2Cluster/conf/Richael_s_pneumo.R")
#   var_order <- c(Main_var,Var_01,Var_02)
#   
#   datesDF <- read_csv(dates_path, col_names = T) %>%
#     dplyr::select(1,any_of(var_order)) 
#   
#   names(datesDF)[1] <- "sampleID"
#   # names(datesDF)[names(datesDF) %in% Var_02]
#   
#   val1 <- which(names(datesDF) %in% Var_02)
#   names(datesDF)[val1] <- "TakenDate"
#   Var_02 <- "TakenDate"
#   
#   source("~/GitHub/SNP2Cluster/EDA_snp_analysis_v2.R")
#   
#   # rm(list=ls())
#   src_path <- file.path("C:/Users/Stanfordk/Documents/GitHub/SNP2Cluster")
#   setwd(src_path)
# }



# # Baby Germs
# 
# Klepp
comparisons <- data.frame(snp=c(11,20,20,25,11),days=c(60,14,60,45,14))
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



# # ACIBA
# comparisons <- data.frame(snp=c(10,14,14,20,25),days=c(14,14,60,60,45))
# for(i in 1:nrow(comparisons)){
#   snpco=comparisons[i,1]
#   daysco=comparisons[i,2]
#   
#   source("~/GitHub/SNP2Cluster/conf/BabyGERMS_kpn_temb.R")
#   
#   source("~/GitHub/SNP2Cluster/EDA_snp_analysis_v2.R")
#   
#   # rm(list=ls())
#   src_path <- file.path("C:/Users/Stanfordk/Documents/GitHub/SNP2Cluster")
#   setwd(src_path)
# }


# # STAAU
# comparisons <- data.frame(snp=c(11,20,20,25,11),days=c(60,14,60,45,14))
# for(i in 1:nrow(comparisons)){
#   snpco=comparisons[i,1]
#   daysco=comparisons[i,2]
#   
#   source("~/GitHub/SNP2Cluster/conf/BabyGERMS_kpn_temb.R")
#   
#   source("~/GitHub/SNP2Cluster/EDA_snp_analysis_v2.R")
#   
#   # rm(list=ls())
#   src_path <- file.path("C:/Users/Stanfordk/Documents/GitHub/SNP2Cluster")
#   setwd(src_path)
# }

# Ecoli
# comparisons <- data.frame(snp=c(11,20,20,25,11),days=c(60,14,60,45,14))
# for(i in 1:nrow(comparisons)){
#   snpco=comparisons[i,1]
#   daysco=comparisons[i,2]
#   
#   source("~/GitHub/SNP2Cluster/conf/BabyGERMS_kpn_temb.R")
#   
#   source("~/GitHub/SNP2Cluster/EDA_snp_analysis_v2.R")
#   
#   # rm(list=ls())
#   src_path <- file.path("C:/Users/Stanfordk/Documents/GitHub/SNP2Cluster")
#   setwd(src_path)
# }



# # FAECAlis
# comparisons <- data.frame(snp=c(11,20,20,25,11),days=c(60,14,60,45,14))
# for(i in 1:nrow(comparisons)){
#   snpco=comparisons[i,1]
#   daysco=comparisons[i,2]
# 
#   source("~/GitHub/SNP2Cluster/conf/BabyGERMS_kpn_temb.R")
# 
#   source("~/GitHub/SNP2Cluster/EDA_snp_analysis_v2.R")
# 
#   # rm(list=ls())
#   src_path <- file.path("C:/Users/Stanfordk/Documents/GitHub/SNP2Cluster")
#   setwd(src_path)
# }

# # FAECium
# comparisons <- data.frame(snp=c(11,20,20,25,11),days=c(60,14,60,45,14))
# for(i in 1:nrow(comparisons)){
#   snpco=comparisons[i,1]
#   daysco=comparisons[i,2]
#   
#   source("~/GitHub/SNP2Cluster/conf/BabyGERMS_kpn_temb.R")
#   
#   source("~/GitHub/SNP2Cluster/EDA_snp_analysis_v2.R")
#   
#   # rm(list=ls())
#   src_path <- file.path("C:/Users/Stanfordk/Documents/GitHub/SNP2Cluster")
#   setwd(src_path)
# }

# Could you do Aciba @Dora Nginza and Staph @Queen Nandi next?