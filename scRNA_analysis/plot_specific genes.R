#small script to visualize some genes of interest

#prepare the environment####
rm(list = ls(all.names = TRUE))
gc()

#load libraries####
library(tidyverse)
source("git_modules/scRNA_tools/scRNA_functions.R")
source("scRNA_analysis/env.R")
source("color_pallets.R")

#get the seurat object
seuratObj <- readRDS("../outputs/scRNA/scran_normalization_no_sct_no_integration/seuratObj_all_timepoints_processed.rds")
#temporary: update colors
seuratObj@misc$colors$atlas.stages.ident <- morphotype_colors

# get the list of meiotic genes from Carolina Catta Preta, 2023
genes_of_interest <- list(
  meiotic_genes_LTR = c(
    HOP1    = "LTRL590_360017200",
    SPO11   = "LTRL590_160011700",
    DMC1    = "LTRL590_350055700",
    HOP2    = "LTRL590_270031700",
    MND1    = "LTRL590_240015900",
    `HAP2-1` = "LTRL590_350009500",
    `HAP2-2` = "LTRL590_360048900"
  )  
)


Ltr_to_Ld <- read_delim(Ltr_to_Ld) %>%
  dplyr::rename(LdCL = source_id, 
                LTR = `Input Ortholog(s)`, 
                ortholog_group = `Ortholog Group`) %>%
  dplyr::select(LdCL, LTR, ortholog_group) %>%
  mutate(LdCL = gsub("_", "-", LdCL)) %>%
  separate_longer_delim(cols = LTR, delim = ",") %>%
  select(LTR, LdCL) %>%
  deframe()

genes_of_interest$meiotic_genes_Ld <-  Ltr_to_Ld[genes_of_interest$meiotic_genes_LTR]

# # get metacyclics markers
metacyclic_markers <- c(
  "sherp",
  "hydrophilic acylated surface protein",
  "gp63",
  "meta "
)

genes_of_interest[["metacyclic-markers"]] <-  seuratObj[["RNA"]]@meta.data %>%
  filter(grepl(paste0(metacyclic_markers, collapse = "|"), tolower(cellranger_name))) %>%
  select(cellranger_name) %>%
  rownames_to_column("rowname") %>%
  select(cellranger_name, rowname) %>%
  deframe()


VlnPlot(seuratObj, features = genes_of_interest$meiotic_genes_Ld, group.by = "seurat_clusters", assay = "RNA")

VlnPlot(seuratObj, features = genes_of_interest$`metacyclic-markers`, group.by = "seurat_clusters", assay = "RNA", alpha = 0.1)