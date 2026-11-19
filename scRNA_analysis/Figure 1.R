#script to plot the figure 1 pannel
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

#set defaults
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("align_plots", "ggalign")
conflict_prefer("unname", "base")
conflict_prefer("crossing", "tidyr")

#get the seurat object
seuratObj <- readRDS("../inputs/scRNA/seurat_objects/sandfly_data/sandfly_processed.rds")

#update colors
seuratObj@misc$colors$assigned_morphotype <- morphotype_colors
seuratObj@misc$colors$pbm_day <- pbm_day_colors

#get cluster markers
cluster_markers <- xlsx::read.xlsx("../outputs/scRNA/sandflies_analysis/Heatmaps/seurat_clusters_markers.xlsx", sheetIndex = 1) %>%
  mutate(feature_id = gsub("_", "-", seurat_id))

#plot cluster umap####
umaps <- list()
umaps[[1]] <- DimPlot(seuratObj, group.by = "pbm_day", cols = seuratObj@misc$colors$pbm_day, label = TRUE) + 
  ggtitle("Clusters") +
  theme(legend.position = "top")

#plot umap per day####
to_plot <- Embeddings(seuratObj, reduction = "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  left_join(
    seuratObj@meta.data %>%
      rownames_to_column("cell")
    , by = "cell")

#create a background dataset for all cells
day_levels <- unique(to_plot$pbm_day)

background <- to_plot %>%
  select(cell, umap_1, umap_2) %>%
  crossing(pbm_day = day_levels) %>%
  mutate(seurat_clusters = NA)

#needed to draw the dashed lines separating the panels
x_rng <- range(to_plot$umap_1)
x_sep <- x_rng[2] + diff(x_rng) * 0.03

#create a dataframe for the labels
labels_df <- to_plot %>%
  group_by(seurat_clusters, pbm_day) %>%
  summarise(umap_1 = median(umap_1), umap_2 = median(umap_2), n_cells = n()) %>%
  group_by(seurat_clusters) %>%
  mutate(prop = n_cells/sum(n_cells)) %>%
  mutate(label = ifelse(prop < 0.1, NA, as.character(seurat_clusters)))

umap_per_day <- ggplot(to_plot, aes(x = umap_1, y = umap_2, color = seurat_clusters))+
  geom_point(data = background, size = 0.00001, show.legend = FALSE, color = "lightgrey") + 
  geom_point(size = 0.00001) + 
  geom_text(data = labels_df, aes(label = label), color = "black", size = 3)+
  geom_vline(xintercept = x_sep, linetype = "dashed",color = "grey40",linewidth = 0.3) +
  facet_wrap(vars(pbm_day), nrow = 1) + 
  scale_color_manual(values = seuratObj@misc$colors$seurat_clusters)+
  guides(name = "Cluster", color = "none")+
  theme_void()+
  theme(legend.spacing.y = grid::unit(0, "cm"), legend.key.height = grid::unit(0.25, "cm"))

#plot clusters per day####
clusters_per_day <- seuratObj@meta.data %>%
  select(pbm_day, seurat_clusters)%>%
  group_by(pbm_day, seurat_clusters)%>%
  summarise(ncells = n(), .groups = "drop") %>%
  complete(pbm_day, seurat_clusters, fill = list(ncells = 0)) %>%
  group_by(pbm_day) %>%
  mutate(percent = ncells/sum(ncells)) %>%
  mutate(percent_text = paste0(round(percent*100, 1), "%"))%>%
  mutate(percent_text = ifelse(percent < 0.1, NA, percent_text)) %>%
  ggplot(aes(x = pbm_day, alluvium = seurat_clusters, stratum = seurat_clusters, y = percent, fill = seurat_clusters, label = percent_text))+
  geom_alluvium(alpha = 0.8, decreasing = FALSE)+geom_text_repel(stat = "stratum", decreasing = FALSE, direction = "y", nudge_x = 0.1)+
  labs(fill = "Cluster", x = NULL, y = "Percentage")+
  scale_fill_manual(values = seuratObj@misc$colors$seurat_clusters)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#plot cluster marker genes####
features_to_plot <- cluster_markers %>%
  mutate(cluster = factor(cluster, levels = levels(seuratObj@meta.data$seurat_clusters))) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  pull(feature_id)

feature_labels <- cluster_markers %>%
  filter(feature_id %in% features_to_plot) %>%
  select(feature_id, cellranger_name) %>%
  deframe()

cluster_markers_plot <- DoHeatmap_plus(
  seuratObj = seuratObj,
  features = features_to_plot,
  group.by = "seurat_clusters",
  group.colors = seuratObj@misc$colors,
  feature_labels = feature_labels,
  feature_label_size = 6,
  feature_labels_max_char = 50
)


#plot specific markers####
features_to_plot <- seuratObj[["RNA"]]@meta.data %>%
  filter(grepl("sherp|dot1", tolower(cellranger_name)))

#format the umap plots
umaps <- lapply(umaps, function(x){
  x <- x +
    scale_x_continuous(breaks = NULL, name = "UMAP 1") +
    scale_y_continuous(breaks = NULL, name = "UMAP 2") +
    ggtitle(NULL)
  return(x)
})


#now build the main figure panel
panel_list <- list()
panel_list[["A"]] <- umaps[[1]]
panel_list[["B"]] <- umap_per_day
panel_list[["C"]] <- clusters_per_day
panel_list[["D"]] <- cluster_markers_plot


#remove label-base alignment
panel_list <- lapply(panel_list, free_lab)

#build sub panels
panel_list[["A-C"]] <- align_plots(panel_list[["A"]], panel_list[["B"]], panel_list[["C"]], nrow = 1, widths = c(0.15, 0.55, 0.2))
panel_list[["A-D"]] <- align_plots(free_align(panel_list[["A-C"]]), free_align(panel_list[["D"]]), ncol = 1, heights = c(0.25, 0.75))


final_plot <- panel_list[["A-D"]] + layout_tags("A")

#save the figure
dir.create("../Figures")
ggsave("../Figures/Figure1.pdf", final_plot, width = 8.27, height = 8.27, scale = 2)
