#this script takes the processed seurat object and map it to the single-cell atlas of Leishmania major.
#clean the environment####
rm(list = ls(all.names = TRUE))
gc()

#get libraries####
library(conflicted)
library(Seurat)
library(tidyverse)
source("git_modules/scRNA_tools/scRNA_functions.R")

#fix package conflicts####
conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")

#get common variables####
source("scRNA_analysis/env.R")
Lmj_to_Lj <- Lmj_to_Lj #here just for readability: it is set by the `env.R` file. 

#set parameters####
Lmj_to_Lj <- xlsx::read.xlsx(Lmj_to_Ld, sheetIndex = 1) #Lmj_to_Ld is set by `env.R`

#get objects#####
Ld_obj <- readRDS("../inputs/scRNA/seurat_objects/sandfly_data/sandfly_processed.rds")
Lmj_obj <- readRDS("../inputs/scRNA/seurat_objects/Lmj_Atlas.rds")

#set metadata to distinct them later####
Ld_obj$species <- "L. donovani"
Lmj_obj$species <- "L. major"

#get anchors####
Anchors <- FindTransferAnchors(query = new_query_obj, reference = ref_obj, normalization.method = normalization_method, reference.reduction = "pca", verbose = FALSE, dims = dims)

