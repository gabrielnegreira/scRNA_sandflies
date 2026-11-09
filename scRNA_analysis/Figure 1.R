#script to plot the figure 1 pannel
#clean the environment
rm(list = ls(all.names = TRUE))
gc()

#libraries#####
library(tidyverse)
library(ggalign)
library(ggalluvial)
library(ggpubr)
library(Seurat)
library(gridExtra)
library(conflicted)
source("git_modules/scRNA_tools/scRNA_functions.R")
source("color_pallets.R")

#set defaults
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("align_plots", "ggalign")

#get the seurat object
seuratObj <- readRDS("../outputs/scRNA/standard_normalization_no_sct_no_integration/seuratObj_all_timepoints_processed.rds")

#temporary: update colors
seuratObj@misc$colors$atlas.stages.ident <- morphotype_colors

#plot umaps####
umaps <- list()
umaps[[1]] <- DimPlot(seuratObj, group.by = "seurat_clusters", cols = seuratObj@misc$colors$seurat_clusters)
umaps[[2]] <- DimPlot(seuratObj, group.by = "atlas.stages.ident", cols = seuratObj@misc$colors$atlas.stages.ident)

#plot the morphotype_to_cluster mapping
cluster_mapping_plot <- seuratObj@meta.data %>%
  group_by(seurat_clusters, atlas.stages.ident) %>%
  summarise(n_cells = n()) %>%
  ggplot(aes(x = seurat_clusters, y = n_cells, fill = atlas.stages.ident))+
  geom_col()+
  scale_fill_manual(values = seuratObj@misc$colors$atlas.stages.ident)

#define a single morphotype for each cluster
morpho_table <- seuratObj@meta.data %>%
  group_by(seurat_clusters, atlas.stages.ident) %>%
  summarise(n_cells = n()) %>%
  group_by(seurat_clusters) %>%
  mutate(n_cells_total = sum(n_cells)) %>%
  mutate(prop = n_cells/n_cells_total) %>%
  filter(prop == max(prop)) 

seuratObj@meta.data <- seuratObj@meta.data %>%
  rownames_to_column("rowname") %>%
  left_join(
    morpho_table %>%
      select(seurat_clusters, atlas.stages.ident) %>%
      rename(assigned_morphotype = atlas.stages.ident)
    , by = "seurat_clusters"
  ) %>%
  column_to_rownames("rowname")

#plot a table with results
header <- tibble(
  row = ".header",
  column = factor(
    c("cluster", "assigned morphotype", "number of cells assigned",
      "cluster total number of cells", "proportion"),
    levels = c("cluster", "assigned morphotype", "number of cells assigned",
               "cluster total number of cells", "proportion")
  ),
  value = c(
    "Cluster", "Assigned morphotype", "Number of cells assigned",
    "Cluster total number of cells", "Proportion"
  ),
  header = TRUE
)

morpho_table_plot <- morpho_table %>%
  mutate(prop = paste0(round(prop * 100, 2), "%")) %>%
  rename(
    cluster = seurat_clusters,
    `assigned morphotype` = atlas.stages.ident,
    `number of cells assigned` = n_cells,
    `cluster total number of cells` = n_cells_total,
    proportion = prop
  ) %>%
  plot_table()
  
#assigned 
umaps[[3]] <- DimPlot(seuratObj, group.by = "assigned_morphotype", cols = seuratObj@misc$colors$atlas.stages.ident)

#format the umap plots
umaps <- lapply(umaps, function(x){
  x <- x +
    scale_x_continuous(breaks = NULL, name = "UMAP 1") +
    scale_y_continuous(breaks = NULL, name = "UMAP 2") +
    ggtitle(NULL)
  
  return(x)
})

#now make the umap showing the assigned morphotypes
plot <- list()
for(Morphotype in levels(seuratObj@meta.data$assigned_morphotype)){
  cells <- seuratObj@meta.data %>%
    filter(assigned_morphotype == Morphotype) %>%
    rownames()
  
  plot[[Morphotype]] <- DimPlot(
    seuratObj, 
    cells.highlight = cells, 
    cols.highlight = seuratObj@misc$colors$atlas.stages.ident[[Morphotype]],
    pt.size = 0.0001
    ) +
    guides(color = "none")+
    ggtitle(Morphotype)+
    theme_void()
}
umap_morphotypes <- align_plots(!!!plot, nrow = 1) + layout_tags(NULL)

#now build the figure panel
plot_list <- list()
plot_list[["A"]] <- ggalign::align_plots(umaps[[1]], umaps[[2]], ncol = 1) + layout_tags(NULL)
plot_list[["B"]] <- cluster_mapping_plot
plot_list[["C"]] <- morpho_table_plot
plot_list[["D"]] <- umap_morphotypes
left <- align_plots(free_align(plot_list[["B"]], axes = "r"), plot_list[["C"]], ncol = 1) 
top <- align_plots(plot_list[["A"]], left, nrow = 1) 
final_plot <- align_plots(top, plot_list[["D"]], ncol = 1, heights = c(0.7, 0.3)) + layout_tags("A")
final_plot

#save the figure
dir.create("../Figures")
ggsave("../Figures/Figure1.pdf", final_plot, width = 8.27, height = 9, scale = 2)
