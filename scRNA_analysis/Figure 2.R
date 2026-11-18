##script to plot the figure 2 pannel
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
conflict_prefer("crossing", "tidyr")

#get the seurat object####
seuratObj <- readRDS("../inputs/scRNA/seurat_objects/sandfly_data/sandfly_processed.rds")

#get cluster markers####
morpho_markers <- xlsx::read.xlsx("../outputs/scRNA/sandflies/Heatmaps/assigned_morphotype_markers.xlsx", sheetIndex = 1) %>%
  mutate(feature_id = gsub("_", "-", seurat_id))

cluster_markers <- xlsx::read.xlsx("../outputs/scRNA/sandflies/Heatmaps/seurat_clusters_markers.xlsx", sheetIndex = 1) %>%
  mutate(feature_id = gsub("_", "-", seurat_id))

#update colors####
seuratObj@misc$colors$assigned_morphotype <- morphotype_colors
seuratObj@misc$colors$pbm_day <- pbm_day_colors

#plot umap with annotations####
morpho_umap <- DimPlot(seuratObj, group.by = "assigned_morphotype", cols = seuratObj@misc$colors$assigned_morphotype)+
  ggtitle(NULL)+
  labs(x = "UMAP 1", y = "UMAP 2")+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)

clusters_umap <- DimPlot(seuratObj, group.by = "seurat_clusters", cols = seuratObj@misc$colors$seurat_clusters, label = TRUE)+
  ggtitle(NULL)+
  labs(x = "UMAP 1", y = "UMAP 2")+
  scale_x_continuous(breaks = NULL)+
  scale_y_continuous(breaks = NULL)+
  guides(color = "none")

#plot the flow between clusters and their morphotype annotation####
lode_colors <- c(seuratObj@misc$colors$assigned_morphotype, seuratObj@misc$colors$seurat_clusters)
annotation_flow_plot <- seuratObj@meta.data %>%
  group_by(seurat_clusters, assigned_morphotype) %>%
  summarise(n_cells = n()) %>%
  ggplot(aes(axis1 = seurat_clusters, axis2 = assigned_morphotype, y = n_cells)) +
  geom_alluvium(aes(fill = seurat_clusters)) +
  geom_stratum(aes(fill = after_stat(stratum))) +
  geom_text(stat = "stratum", aes(label = after_stat(str_remove(stratum, "\\.\\d+$")))) +
  scale_fill_manual(values = lode_colors) +
  labs(x = NULL, y = "Number of cells") +
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#temporary: plot heatmap of marker genes for annotated clusters
#plot cluster marker genes####
features_to_plot <- morpho_markers %>%
  mutate(cluster = factor(cluster, levels = levels(seuratObj@meta.data$assigned_morphotype))) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  pull(feature_id)

features_to_plot <- cluster_markers %>%
  mutate(cluster = factor(cluster, levels = levels(seuratObj@meta.data$seurat_clusters))) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  pull(feature_id)

feature_labels <- cluster_markers %>%
  filter(feature_id %in% features_to_plot) %>%
  select(feature_id, cellranger_name) %>%
  deframe()

markers_plot <- DoHeatmap_plus(
  seuratObj = seuratObj,
  features = features_to_plot,
  group.by = c("seurat_clusters", "assigned_morphotype"),
  group.colors = seuratObj@misc$colors,
  feature_labels = feature_labels,
  feature_label_size = 9
)

#monut the pannel####
panel_list <- list()
panel_list[["A"]] <- clusters_umap
panel_list[["B"]] <- morpho_umap
panel_list[["C"]] <- annotation_flow_plot
panel_list[["D"]] <- markers_plot

#remove label-base alignment####
panel_list <- lapply(panel_list, free_lab)

#build the panel ####
panel_list[["A-C"]] <- align_plots(panel_list[["A"]], panel_list[["B"]], panel_list[["C"]], nrow = 1)
panel_list[["A-D"]] <- align_plots(free_align(panel_list[["A-C"]]), free_align(panel_list[["D"]]), ncol = 1, heights = c(0.33, 0.66))

#add tags ####
final_plot <- panel_list[["A-D"]] + layout_tags("A")

#export it ####
dir.create("../Figures")
ggsave("../Figures/Figure2.pdf", final_plot, width = 8.27, height = 8.27, scale = 2)



