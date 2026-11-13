#temporary script to test the `DoHeatmap_plus` function
source("git_modules/scRNA_tools/scRNA_functions.R")
library(tidyverse)

seuratObj <- readRDS("../inputs/scRNA/seurat_objects/seuratObj_all_timepoints_processed.rds")
marker_genes <- xlsx::read.xlsx("../outputs/scRNA/Heatmaps/assigned_morphotype_final_markers.xlsx", sheetIndex = 1) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  mutate(id = gsub("_", "-", cellranger_id))

genes_to_plot <- marker_genes$id
gene_names <- marker_genes %>%
  ungroup() %>%
  dplyr::select(id, cellranger_name) %>%
  deframe()

#set the variables
cells = NULL
group.by = "ident"
layer = "scale.data"
assay = NULL
show_feature_labels = NULL
feature_labels = NULL
feature_labels_max_char = 40
feature_label_size = 11
group.colors = NULL
hclust_method = "complete"
expression_colors = NULL
raster = NULL
raster_resolution = 300
seuratObj <- seuratObj
features = genes_to_plot
feature_labels = gene_names
group.by = c("assigned_morphotype_final", "seurat_clusters")
group.colors = seuratObj@misc[["colors"]]

DoHeatmap_plus(
  seuratObj, 
  features = genes_to_plot, 
  feature_labels = gene_names, 
  group.by = c(
    #"assigned_morphotype_final"#, 
    "seurat_clusters"
    ),
  group.colors = seuratObj@misc$colors
  )