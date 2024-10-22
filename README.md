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
  1. Pairwise single nucleotide polymorphism (SNP) distance matrix 
  2. SNP cut-off and interval between specimen collection dates (defaults: 20 SNPs and 14 days), 
  3. Multi-locus sequence type (MLST) profiles
  4. Epidemiological data

Closely-related isolates are initially grouped together in clusters based on an enhanced K-means clustering method that employs a silhouette score and 500 bootstraps to determine the optimal K (which defines the maximum number of clusters) **K-means clusters**. Custom R functions are applied to each pre-grouped cluster to generate SNP cluster chains based on the provided SNP cut-off to make a SNP cluster, and resetting when SNP threshold is exceeded to make the next cluster. Each isolate is only assigned once to a cluster. 

If MLST profiles and epidemiological data are provided, the final transmission clusters are generated in the context of sequence types and epidemiological timeline. Publication-ready graphs are automatically generated for visualization of the integrated SNP/Epi transmission clusters, including a heat-map, minimum-spanning tree and scatter plots.
