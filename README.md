[![DOI](https://zenodo.org/badge/807615693.svg)](https://doi.org/10.5281/zenodo.13976252)

# SNP2Cluster

> ### Core SNP-based clustering for enhanced transmission cluster detection in outbreak scenarios

SNP2Cluster is a K-means clustering-based method that identifies transmission clusters by integrating genomic and epidemiological data. 

Four types of cluster analyses are possible using SNP2Cluster:

1. Core SNP cluster analysis - using SNP data only
2. Transmission cluster analysis which integrates genomic and epi data:
    1. Facility level analysis - STs and clusters visualized strictly per facility
    1. Area/location level analysis - STs and facilities grouped by location
    1. Community level analysis - facilities grouped by ST
  
## Methods overview
To perform enhanced transmission cluster analysis, paths to the following information/data sources should be provided in a configuration file: 
  1. Epidemiological data file including collection dates and facility 
     information
  2. Pairwise single nucleotide polymorphism (SNP) distance matrix 
  3. Multi-locus sequence type (MLST) profiles
  4. Output directory

Closely-related isolates are initially grouped together in clusters based on an enhanced K-means clustering method that employs a silhouette score and 500 bootstraps to determine the optimal K (which defines the maximum number of clusters) **K-means clusters**. Custom R functions are applied to each pre-grouped cluster to generate SNP cluster chains based on the provided SNP cut-off to make a SNP cluster, and resetting when SNP threshold is exceeded to make the next cluster. Each isolate is only assigned once to a cluster **Core SNP cluster analysis**. 

If MLST profiles and epidemiological data are provided, the final transmission clusters are generated in the context of sequence types and epidemiological timeline **Transmission clusters**. Publication-ready graphs are automatically generated for visualization of the integrated SNP/Epi transmission clusters, including a heat-map, minimum-spanning tree and scatter plots.

## Parameter setting and usage
**1. Provide paths to the input files in the configuration file - refer to the config_file_template in the conf folder**

> Refer to the **example-data** folder to see examples of input files

```
# Set paths to files ------------------------------------------------------

dates_path = "./example-data/example_metadata.csv"     # Epidemiological data file including collection dates and facility information
filepath = "./example-data/coreSNPmatrix.co.csv"       # Pairwise single nucleotide polymorphism (SNP) distance matrix
mlst_profile = "./example-data/05.mlst.xlsx"           # Multi-locus sequence type (MLST) profiles
out_dir <- "./example-output"                          # Output directory

```

**2. The following variables should be defined in the configuration file as well**

```
# Define variables --------------------------------------------------------

# The first column of the metadata file (dates_path) should have the sample_ids
# The other variables such as Facility, collection dates etc should be assigned to the variables below:

Main_var = ""   # Mandatory main variable e.g. Hospital or Facility etc.
Var_01 =  ""    # Optional second variable e.g. Ward_name, Ward_type etc. 
Var_02 = ""     # Mandatory variable for specimen collection dates
clust_type = "" # "Core" or "Transmission"
snpco = 20      # set SNP cut-off     (Set to 20 by default)

if(clust_type == "Transmission"){
  trans_lvl = "Facility" #"Community" # Community or "Facility" #Default 
                         # facility
  daysco = 45 # Time interval in days (Default: 45 days)
}

```
**3. Specify format of collection dates in the epi data file**

```
# Specify format of collection dates --------------------------------------

lubri_fmt <- "ymd" 
# collection date format
# options include: dym, dmy, ymd, ydm, etc.. -- based on the lubridate package

```

>  _Save the configuration file in the **conf folder**_


**4. Set additional parameters in the execution file**

_**a. Override the defaults for snp threshold and days interval and provide your preferred in intervals**_
```
rm(list=ls())
snpco=20    # Preferred SNP threshold
daysco=45   # Time interval in days

# Optionally a vector of SNP threshols can be provided together with matching
# day intervals in a second vector if you want to perform multiple comparisons

# snpco=c(11,20,20,25,11)
# snpco=c(60,14,60,45,14)

```

_**b. Provide name of the configuration file saved in the conf folder**_

```
conf_file="config_file_template.R"  # config file should be saved in the conf folder    

```

**5. Run the analysis**
```
source("run_analysis.R")
```

## Reference (Citation)
Kwenda, S., Shuping, L., Mashau, R., Ismail, H., & Govender, N. P. (2024). SNP2Cluster: A core SNP and K-means clustering-based tool for enhanced transmission cluster detection in outbreak scenarios (v0.5.3). Klebsiella Epidemiology and Biology Symposium 2024 (KLEBS2024), Institut Pasteur, Paris, France. Zenodo. [https://zenodo.org/records/14060296](https://zenodo.org/records/14060296)

Shuping, L., Ismail, H., Mashau, R., Kwenda, S., Holt, K. E., Magobo, R. E., Perovic, O., Meiring, S. T., Quan, V. C., & Govender, N. P. (2025). Enhanced detection of neonatal invasive infection clusters in South Africa using epidemiological and genomic surveillance data. medRxiv. https://doi.org/10.1101/2025.11.10.25339895
