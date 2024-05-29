# library(devtools)
# library(bitrugs)
library(plyr)
library(tidyverse)
library(readr)
library(adegenet)
library(outbreaker2)
library(TransPhylo)
library(ape)
library(igraph)
library(lubridate)

### outbreaker2 examples ####
library(ape)
library(outbreaker2)
str(fake_outbreak)
fake_outbreak$dna
fake_outbreak$onset
fake_outbreak$sample
fake_outbreak$ances
fake_outbreak$ctd
## adegenet
dat <- haploGen(seq.l=1e4, repro=function(){sample(1:4,1)}, gen.time=1, t.max=3)
dat
str(dat)
dat <- read.csv(system.file("files/pdH1N1-data.csv",package="adegenet"))
xy <- cbind(dat$lon, dat$lat)
temp <- as.matrix(dist(xy))
M <- 1* (temp < 1e-10)
## SEQTRACK ANALYSIS
res <- seqTrack(dat, mu=0.0001, haplo.length=1e4) 
dat$date <- as.POSIXct(dat$date)
#############################


# devtools::install_github('xavierdidelot/TransPhylo')
help(TransPhylo)
vignette("TransPhylo")

# Set paths ---------------------------------------------------------------
path="/Users/stanfordk/Documents/Data-Delivery/NICD/SCF/CHARM/Rindidzani/Results-CRE-samples-45/iqtree-rpt"
aln_path = "/Users/stanfordk/Documents/Data-Delivery/NICD/SCF/CHARM/Rindidzani/snippy/snippy_edt_aln.fasta"
dates_path = "/Users/stanfordk/Documents/Data-Delivery/NICD/SCF/CHARM/Rindidzani/snippy/date_file.tsv"
list.files(path = path)

coll_dates <- read_tsv(dates_path, col_names = F)

# Read in DNA alignment sequences -----------------------------------------
aln_seq <- read.dna(aln_path, format = "fasta")
aln_names <- labels(aln_seq)
# labels(aln_seq) <- str_remove_all(labels(aln_seq),".fasta")
aln_dist <- dist.dna(aln_seq)
str(aln_dist)



# Read tree file ----------------------------------------------------------
##
phy<-read.tree(text='((4:1.257652937,(1:1.231048819,5:1.519248672):0.303038892):0.784065883,(3:1.643413444,(6:0.656820028,2:0.007344035611):0.7562780805):1.293120815);')
plot(phy)
axisPhylo(backward = F)

t1 <- ptreeFromPhylo(ape::rtree(5),2020)
plot(t1)
##
phy <- read.nexus(file.path(path,"snippy_edt.aln.timetree.nex"))
plot(phy)
axisPhylo(backward = F)

ptree<-ptreeFromPhylo(phy,dateLastSample=2021)
plot(ptree)

####################
