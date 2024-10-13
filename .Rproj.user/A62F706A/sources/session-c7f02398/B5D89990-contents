# ###################################### NOTES ################################################

# The first 3 columns of the metadata file (dates_path) should have the following - 
# in the order specified below:
#   1. sample_id
#   2. collection_date
#   3. facility_name/hospital_name/community codes/regions/metros
#   4. Ward_name (or ward type) if 3 is hospital or any other variable

################################### SET PARAMETERS HERE ########################################

# NOTES -------------------------------------------------------------------

# The first 3 columns of the metadata file (dates_path) should have the following - 
# in the order specified below:
#   1. sample_id
#   2. collection_date
#   3. facility_name/hospital_name/community codes/regions/metros
#   4. Ward_name (or ward type) if 3 is hospital or any other variable


# Set paths to files ------------------------------------------------------

# # # Richael - S pneumoniae
dates_path = "D:/Terra-Informatix/Richael-Ghana/Spneumo-WGS/Genome_Pneumo_Metadata.csv"    #Metatdata-02_STs.csv
filepath = "D:/Terra-Informatix/Richael-Ghana/Spneumo-WGS/snps.coresites.01.SNPmatrix.csv"
moltenpath = "D:/Terra-Informatix/Richael-Ghana/Spneumo-WGS/snps.coresites.01.SNPmatrixmolten.csv"
mlst_profile = "D:/Terra-Informatix/Richael-Ghana/Spneumo-WGS/sequence_types.csv"
out_dir <- "D:/Terra-Informatix/Richael-Ghana/Spneumo-WGS/clusters-05"

# Define variables --------------------------------------------------------

Main_var = "Location"
# Possible values for Main_var based on provided metadata
# Community, Hospital, Facility

Var_01 = "Facility"
# Possible values for Var_01 can be:
# Ward, WardType, Facilities within a community (if Main_var == "Community")


Var_02 = "SamplingDate"
# Specify name of column with collection dates

clust_type = "Transmission" # Core or Transmission

if(clust_type == "Transmission"){
  # Type of transmission dynamics (Community vs Facility vs Hospital)
  # transmission_type = "facility" # Facility level analysis will be performed per facility. Facility can be a school, hopsital etc
  # transmission_type = "facility" # Hospital level analysis will be performed per hospital. Ward information should be provided
  transmission_type = Main_var #"community" # Community level analysis will be performed by ST
  #snpco=11  # set SNP cut-off     (Set to 20 by default)
  #daysco=14 # Set number of days (set to 14 by default)
  
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





# # Type of transmission dynamics (Community vs Facility vs Hospital)
# transmission_type = "facility" # Facility level analysis will be performed per facility. Facility can be a school, hopsital etc
# # transmission_type = "facility" # Hospital level analysis will be performed per hospital. Ward information should be provided
# # transmission_type = "community" # Community level analysis will be performed by ST
# snpco=20  # set SNP cut-off     (Set to 20 by default)
# daysco=21 # Set number of days (set to 14 by default)
# 
# # # Transmission network input files
# # path="D:/Terra-Informatix/Richael-Ghana/Spneumo-WGS/transmission-analysis"
# # aln_path = "D:/Terra-Informatix/Richael-Ghana/Spneumo-WGS/snps.coresites.01.snp_sites.aln"
# # snp_dist <- moltenpath
# # core_snps = filepath
# 
# 
# # collection date format
# lubri_fmt <- "mdy" # options include: dym, dmy, ymd, ydm, etc.. -- based on the lubridate package
# 
# refDate="1800-01-01"
# refST=NA
# # NOTES
# 
# # The first 3 columns of the metadata file (dates_path) should have the following - 
# # in the order specified below:
# #   1. sample_id
# #   2. collection_date
# #   3. facility_name/hospital_name/community codes/regions/metros
# #   4. Ward_name (or ward type) if 3 is hospital or any other variable
# 
# # 
# # 
# # 
# # 
# #




