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

# Baby GERMS Study --------------------------------------------------------------

# KLEBS paths -------------------------------------------------------------------
# # TEMBISA
# dates_path = "E:/projects/Baby-Germs/KLEPP/KLEPP/TEMBISA/transmission/TEMBISA_ward.csv"
# filepath = "E:/projects/Baby-Germs/KLEPP/KLEPP/TEMBISA/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/KLEPP/KLEPP/TEMBISA/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/KLEPP/KLEPP/TEMBISA/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/KLEPP/TEMBISA/clusters2409"

# # DORA
# dates_path = "E:/projects/Baby-Germs/KLEPP/KLEPP/DORA-NGINZA/transmission/DORA_ward.csv"
# filepath = "E:/projects/Baby-Germs/KLEPP/KLEPP/DORA-NGINZA/DNH/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/KLEPP/KLEPP/DORA-NGINZA/DNH/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/KLEPP/KLEPP/DORA-NGINZA/DNH/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/KLEPP/DORA-NGINZA/clusters2409"

# # hospital="KLERKSDORP"
# dates_path = "E:/projects/Baby-Germs/KLEPP/KLEPP/KLERKSDORP/transmission/KLERKSDORP_ward.csv"
# filepath = "E:/projects/Baby-Germs/KLEPP/KLEPP/KLERKSDORP/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/KLEPP/KLEPP/KLERKSDORP/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/KLEPP/KLEPP/KLERKSDORP/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/KLEPP/KLERKSDORP/clusters2409"

# # MANKWENG
# dates_path = "E:/projects/Baby-Germs/KLEPP/KLEPP/MANKWENG/transmission/MANKWENG_ward.csv"
# filepath = "E:/projects/Baby-Germs/KLEPP/KLEPP/MANKWENG/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/KLEPP/KLEPP/MANKWENG/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/KLEPP/KLEPP/MANKWENG/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/KLEPP/MANKWENG/clusters2409"

# # # QUEEN-NANDI-REGIONAL
# dates_path = "E:/projects/Baby-Germs/KLEPP/KLEPP/QUEEN-NANDI-REGIONAL/transmission/QUEEN_ward.csv"
# filepath = "E:/projects/Baby-Germs/KLEPP/KLEPP/QUEEN-NANDI-REGIONAL/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/KLEPP/KLEPP/QUEEN-NANDI-REGIONAL/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/KLEPP/KLEPP/QUEEN-NANDI-REGIONAL/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/KLEPP/QUEEN-NANDI-REGIONAL/clusters2409"

# # ROB-FERREIRA
# dates_path = "E:/projects/Baby-Germs/KLEPP/KLEPP/ROB-FERREIRA/transmission/ROB_ward.csv"
# filepath = "E:/projects/Baby-Germs/KLEPP/KLEPP/ROB-FERREIRA/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/KLEPP/KLEPP/ROB-FERREIRA/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/KLEPP/KLEPP/ROB-FERREIRA/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/KLEPP/ROB-FERREIRA/clusters2409"


# ACIBA paths -------------------------------------------------------------
# # # DORA
# dates_path = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/DORA/transmission/DORA_ward.csv"
# filepath = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/DORA/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/DORA/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/DORA/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/ACIBA/DORA-NGINZA/clusters2409"

# # # # KLERKSDORP
# dates_path = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/KLERKSDORP/transmission/KLERKSDORP_ward.csv"
# filepath = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/KLERKSDORP/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/KLERKSDORP/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/KLERKSDORP/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/ACIBA/KLERKSDORP/clusters2409"

# # # MANKWENG
# dates_path = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/MANKWENG/transmission/MANKWENG_ward.csv"
# filepath = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/MANKWENG/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/MANKWENG/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/MANKWENG/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/ACIBA/MANKWENG/clusters2409"

# # QUEEN
# dates_path = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/QUEEN/transmission/QUEEN_ward.csv"
# filepath = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/QUEEN/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/QUEEN/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/QUEEN/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/ACIBA/QUEEN-NANDI-REGIONAL/clusters2409"

# # ROB
# dates_path = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/ROB/transmission/ROB_ward.csv"
# filepath = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/ROB/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/ROB/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/ROB/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/ACIBA/ROB-FERREIRA/clusters2409"
 
# # TEMBISA
# dates_path = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/TEMBISA/transmission/TEMBISA_ward.csv"
# filepath = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/TEMBISA/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/TEMBISA/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/ACIBA/downstream-analysis/ACIBA/TEMBISA/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/ACIBA/TEMBISA/clusters2409"



# STAPH -------------------------------------------------------------------
# # DORA
# dates_path = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/DORA/transmission/DORA_ward.csv"
# filepath = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/DORA/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/DORA/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/DORA/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/STAAU/DORA-NGINZA/clusters2409"

# # # # KLERKSDORP
# dates_path = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/KLERKSDORP/transmission/KLERKSDORP_ward.csv"
# filepath = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/KLERKSDORP/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/KLERKSDORP/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/KLERKSDORP/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/STAAU/KLERKSDORP/clusters2409"

# # # MANKWENG
# dates_path = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/MANKWENG/transmission/MANKWENG_ward.csv"
# filepath = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/MANKWENG/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/MANKWENG/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/MANKWENG/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/STAAU/MANKWENG/clusters2409"

# # QUEEN
# dates_path = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/QUEEN/transmission/QUEEN_ward.csv"
# filepath = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/QUEEN/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/QUEEN/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/QUEEN/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/STAAU/QUEEN-NANDI-REGIONAL/clusters2409"

# # ROB
# dates_path = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/ROB/transmission/ROB_ward.csv"
# filepath = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/ROB/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/ROB/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/ROB/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/STAAU/ROB-FERREIRA/clusters2409"

# # TEMBISA
# dates_path = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/TEMBISA/transmission/TEMBISA_ward.csv"
# filepath = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/TEMBISA/transmission/coreSNPmatrix.co.csv"
# moltenpath = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/TEMBISA/transmission/coreSNPmatrixmolten.csv"
# mlst_profile = "E:/projects/Baby-Germs/STAAU/downstream-analysis/STAAU/TEMBISA/transmission/05.mlst.xlsx"
# out_dir <- "D:/Terra-Informatix/baby-germs/STAAU/TEMBISA/clusters2409"







# Define variables --------------------------------------------------------

Main_var = "Hospital"
# Possible values for Main_var based on provided metadata
# Community, Hospital, Facility

Var_01 = "WardType"
# Possible values for Var_01 can be:
# Ward, WardType, Facilities within a community (if Main_var == "Community")


Var_02 = "TakenDate"
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

lubri_fmt <- "ymd" 
# collection date format
# options include: dym, dmy, ymd, ydm, etc.. -- based on the lubridate package

refDate="1800-01-01"
refST=NA











# # collection date format
# lubri_fmt <- "ymd" # options include: dym, dmy, ymd, ydm, etc.. -- based on the lubridate package
# 

# # # # 1. Sensitivity analysis for Klebsiella
# # # # ·   SNP <11 + 60 days
# # # # ·   SNP <20 + 14 days
# # # # ·   SNP <20 + 60 days
# # # # ·   SNP and days according to Kats paper
# # # # .   SNP <=25 + 45 days 
# # 

# # # # # # KLEP set ST for reference
# refST=11
# refDate="2011-01-01"

# # # # KLEP Get different combinations to use for SNP and EPI days
# comparisons <- data.frame(snp=c(11,20,20,25,11),days=c(60,14,60,45,14))

