# SNP2Cluster

> ### Core SNP-based clustering for enhanced transmission cluster detection in outbreak scenarios

SNP2Cluster is a K-means clustering-based method that identifies transmission clusters by integrating genomic and epidemiological data. Four types of cluster analyses are possible using SNP2Cluster:

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
**1. The following variables should be defined in the configuration file as well**

```
# Define variables --------------------------------------------------------

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
**2. Specify format of collection dates in the epi data file**

  _Save the configuration file in the **conf folder**_
```
# Specify format of collection dates --------------------------------------

lubri_fmt <- "ymd" 
# collection date format
# options include: dym, dmy, ymd, ydm, etc.. -- based on the lubridate package

```
**3. Set additional parameters in the execution file**

```

```

**4. Run the analysis**
```
source("run_analysis.R")
```
