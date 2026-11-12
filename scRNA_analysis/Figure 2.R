# This script plots Figure 2 of the paper
#clean the environment
rm(list = ls(all.names = TRUE))
gc()

#libraries#####
library(tidyverse)
library(ggalign)
library(ggrepel)
library(ggalluvial)
library(ggpubr)
library(Seurat)
library(stringr)
library(conflicted)
source("git_modules/scRNA_tools/scRNA_functions.R")
source("color_pallets.R")

#set defaults####
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("align_plots", "ggalign")
conflict_prefer("unname", "base")

#get the seurat object####
seuratObj <- readRDS("../inputs/scRNA/seurat_objects/seuratObj_all_timepoints_processed.rds")


#

