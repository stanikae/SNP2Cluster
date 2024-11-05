################################### SET PARAMETERS HERE ########################################

# NOTES -------------------------------------------------------------------

# The first column of the metadata file (dates_path) should have the sample_ids
# The other variables such as Facility, collection dates etc should be assigned to the variables below:

# Main_var # Mandatory: The main variable e.g. Facility should be assigned to this variable
# Var_01 = # Optional: If e.g main variable is for hospitals, this variable can be for Wards, Ward type etc.
# Var_02 = # Madatory for Transmission analysis - collection dates should be assigned to this variable. 
#          # Not required for Core SNP cluster analysis


# The variables below are used to specify the type of cluster analysis to be performed, i.e. Core vs Transmission clusters
# 1. Core clusters are only calculated based on SNP distances and require a SNP threshold to be specified. 
#    By default, SNP threshold will be set to 20.
# 2. Transmission clusters are calculated using both SNP thresholds and epi data such as facility/location, time intervals, collection dates
#    and are visualized in the context of sequence types (ST)

# clust_type = "Transmission" # Core or Transmission
# snpco = # SNP threshold
# daysco = # Time interval in days (Default: 45 days)

# The variable below is used to specify if the Transmission analysis is at the level of the facility or community,
# if community is set, then the visualization of the transmission clusters will be by ST. 
# Community level analysis looks at transmission across multiple facilities within a specified location or area in the 
# context of STs and in a given time period

# trans_lvl = "Facility" #"Community" # Community or "Facility" #Default facility


# Lastly you have to specify date format of collection date based on the lubridate package's date formats:
# lubri_fmt = # e.g. "mdy" for Month-Day-Year date format (10-24-2024) or "ymd" for Year-Month-Day (2024-10-24), etc.



# Set paths to files ------------------------------------------------------

dates_path = "./example-data/example_metadata.csv"
filepath = "./example-data/coreSNPmatrix.co.csv"
mlst_profile = "./example-data/05.mlst.xlsx"
out_dir <- "./example-output"

out_dir <- normalizePath(out_dir)
# Define variables --------------------------------------------------------

Main_var = "FacilityName"
Var_01 =  "WardType"
Var_02 = "TakenDate"
clust_type = "Transmission" # Core or Transmission
snpco = 20 # set SNP cut-off     (Set to 20 by default)

if(clust_type == "Transmission"){
  trans_lvl = "Facility" #"Community" # Community or "Facility" #Default facility
  daysco = 45 # Time interval in days (Default: 45 days)
}


# Specify format of collection dates --------------------------------------

lubri_fmt <- "ymd" 
# collection date format
# options include: dym, dmy, ymd, ydm, etc.. -- based on the lubridate package

refDate="1800-01-01"
refST=NA
