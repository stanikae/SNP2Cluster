
# NOTES -------------------------------------------------------------------

# The first 3 columns of the metadata file (dates_path) should have the following - 
# in the order specified below:
#   1. sample_id
#   2. collection_date
#   3. facility_name/hospital_name/community codes/regions/metros
#   4. Ward_name (or ward type) if 3 is hospital or any other variable



# Set paths to files ------------------------------------------------------


# # Thabo-UP P. aeroginosa
dates_path = "E:/projects/Thabo-ST242/ST242-20230904/metadata_file2.xlsx"
filepath = "E:/projects/Thabo-ST242/ST242-20230904/snpdist/coreSNPmatrix.co.csv"
moltenpath = "E:/projects/Thabo-ST242/ST242-20230904/snpdist/coreSNPmatrixmolten.csv"
mlst_profile = "E:/projects/Thabo-ST242/ST242-20230904/sequence_types.xlsx"
out_dir <- "E:/projects/Thabo-ST242/ST242-20230904/cluster-analysis/clusters"
# # Transmission network input files
# path="E:/projects/Thabo-ST242/ST242-20230904/cluster-analysis/minimum-spanning-tree"
# aln_path = "E:/projects/Thabo-ST242/ST242-20230904/snpsites/snps.co.snp_sites.aln"
# snp_dist <- moltenpath
# core_snps = filepath


# Define variables --------------------------------------------------------

Main_var = "ST"
# Possible values for Main_var based on provided metadata
  # Community, Hospital, Facility

Var_01 = "Country"
# Possible values for Var_01 can be:
# Ward, WardType, Facilities within a community (if Main_var == "Community")


Var_02 = "collection_year"
# Specify name of column with collection dates

clust_type = "Core" # Core or Transmission

if(clust_type == "Transmission"){
  # Type of transmission dynamics (Community vs Facility vs Hospital)
  # transmission_type = "facility" # Facility level analysis will be performed per facility. Facility can be a school, hopsital etc
  # transmission_type = "facility" # Hospital level analysis will be performed per hospital. Ward information should be provided
  transmission_type = Main_var #"community" # Community level analysis will be performed by ST
  snpco=20  # set SNP cut-off     (Set to 20 by default)
  daysco=60 # Set number of days (set to 14 by default)
  
}else{
  transmission_type = Main_var # Value for Main_var
  snpco=20  # set SNP cut-off     (Set to 20 by default)
}




# Specify format of collection dates --------------------------------------

lubri_fmt <- "mdy" 
# collection date format
# options include: dym, dmy, ymd, ydm, etc.. -- based on the lubridate package

refDate="1800-01-01"
refST=NA




# 
# 
# 
# 
# 





