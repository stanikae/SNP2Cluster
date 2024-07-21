# ###################################### NOTES ################################################

# The first 3 columns of the metadata file (dates_path) should have the following - 
# in the order specified below:
#   1. sample_id
#   2. collection_date
#   3. facility_name/hospital_name/community codes/regions/metros
#   4. Ward_name (or ward type) if 3 is hospital or any other variable

################################### SET PARAMETERS HERE ########################################

# BabyGERMS
# # TEMBISA
dates_path = "E:/projects/Baby-Germs/KLEPP/KLEPP/TEMBISA/transmission/TEMBISA_ward.csv"
filepath = "E:/projects/Baby-Germs/KLEPP/KLEPP/TEMBISA/transmission/coreSNPmatrix.co.csv"
moltenpath = "E:/projects/Baby-Germs/KLEPP/KLEPP/TEMBISA/transmission/coreSNPmatrixmolten.csv"
mlst_profile = "E:/projects/Baby-Germs/KLEPP/KLEPP/TEMBISA/transmission/05.mlst.xlsx"
out_dir <- "D:/Terra-Informatix/baby-germs/KLEPP/TEMBISA/clusters"


# Type of transmission dynamics (Community vs Facility vs Hospital)
# transmission_type = "facility" # Facility level analysis will be performed per facility. Facility can be a school, hopsital etc
transmission_type = "hospital" # Hospital level analysis will be performed per hospital. Ward information should be provided
# transmission_type = "community" # Community level analysis will be performed by ST
snpco=25#11  # set SNP cut-off     (Set to 20 by default)
daysco=45#60 #14 # Set number of days (set to 14 by default)


# # Transmission network input files
# path="E:/projects/Baby-Germs/KLEPP/KLEPP/TEMBISA/cluster-analysisZ01/transmission-analysis"
# aln_path = "E:/projects/Baby-Germs/KLEPP/KLEPP/TEMBISA/transmission/snps.co.snp_sites.aln"
# snp_dist <- moltenpath
# core_snps = "E:/projects/Baby-Germs/KLEPP/KLEPP/TEMBISA/transmission/coreSNPmatrix.co.csv"

# collection date format
lubri_fmt <- "ymd" # options include: dym, dmy, ymd, ydm, etc.. -- based on the lubridate package


# # # # 1. Sensitivity analysis for Klebsiella
# # # # ·   SNP <11 + 60 days
# # # # ·   SNP <20 + 14 days
# # # # ·   SNP <20 + 60 days
# # # # ·   SNP and days according to Kats paper
# # # # .   SNP <=25 + 45 days 
# # 

# # # # # KLEP set ST for reference
refST=11
refDate="2011-01-01"

# # # # KLEP Get different combinations to use for SNP and EPI days
# comparisons <- data.frame(snp=c(11,20,20,25,11),days=c(60,14,60,45,14))

