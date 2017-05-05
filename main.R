# Setup : Libraries
library(plyr)
library(dplyr)
library(lattice)
library(ggplot2)
library(xtable)
library(gridExtra)
library(stringr)
library(topGO)
library(ALL)
library(ssh.utils)

# source("http://bioconductor.org/biocLite.R") # In case bioMart packages aren't installed
# biocLite("topGO")

# Setup Environment
setwd("/home/bcallaghan/Projects/Onenineteen/")

# Source scripts
source("func.R")
source("0-func.r")
source("load.R")
source("clean.R")
source("do.R")
