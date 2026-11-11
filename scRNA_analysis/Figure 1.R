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

#get the seurat object
seuratObj <- readRDS("../inputs/scRNA/seurat_objects/seuratObj_all_timepoints_processed.rds")

#temporary: test different cluster resolutions
#seuratObj$seurat_clusters <- seuratObj$RNA_snn_res.0.5
#seuratObj@misc$colors$seurat_clusters <- create_colors(levels(seuratObj$seurat_clusters))

#plot umaps####
umaps <- list()
umaps[[1]] <- DimPlot(seuratObj, group.by = "seurat_clusters", cols = seuratObj@misc$colors$seurat_clusters, label = TRUE) + ggtitle("Clusters")
umaps[[2]] <- DimPlot(seuratObj, group.by = "atlas.stages.ident", cols = seuratObj@misc$colors$atlas.stages.ident) + ggtitle("Raw Mapping to Atlas")

#plot the morphotype_to_cluster mapping
cluster_mapping_plot <- seuratObj@meta.data %>%
  group_by(seurat_clusters, atlas.stages.ident) %>%
  summarise(n_cells = n()) %>%
  ggplot(aes(x = seurat_clusters, y = n_cells, fill = atlas.stages.ident))+
  geom_col()+
  scale_fill_manual(values = seuratObj@misc$colors$atlas.stages.ident)

#plot umap with final assigned morphotypes
umaps[[3]] <- DimPlot(seuratObj, group.by = "assigned_morphotype_final", cols = seuratObj@misc$colors$assigned_morphotype_final) + ggtitle("final maping")

#plot the flow for anotation####
##set the columns to plot
flow_columns <- c("atlas.stages.ident","seurat_clusters", "subcluster", "assigned_morphotype_final")

##set subcluster factor levels
subcluster_levels <- seuratObj@meta.data %>%
  distinct(subcluster) %>%
  mutate(
    cluster_id = as.integer(str_extract(subcluster, "^\\d+")),
    subcluster_id = as.integer(str_extract(subcluster, "(?<=_)\\d+$")),
    subcluster_id = coalesce(subcluster_id, 0L)
  ) %>%
  arrange(cluster_id, subcluster_id) %>%
  pull(subcluster)

##set color vector for strata
flow_colors <- seuratObj@misc$colors[flow_columns]

###use original seurat clusters for all subclusters
for(Cluster in names(flow_colors$subcluster)){
  original_cluster <- gsub("_.*", "", Cluster)
  flow_colors$subcluster[Cluster] <- flow_colors$seurat_clusters[original_cluster]
}

###merge all colors into a single vector
flow_colors_vector <- do.call(c, unname(flow_colors))
flow_colors_vector <- flow_colors_vector[!duplicated(names(flow_colors_vector))]

##create a long-format dataframe for the alluvial plot
lodes <- seuratObj@meta.data %>%
  mutate(subcluster = factor(subcluster, levels = subcluster_levels)) %>%
  group_by(across(all_of(flow_columns))) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  to_lodes_form(axes = seq_along(flow_columns), discern = TRUE) %>%
  mutate(stratum_label = str_remove(as.character(stratum), "\\.\\d+$"))


##update the color vector
lode_colors <- setNames(flow_colors_vector[lodes$stratum_label], as.character(lodes$stratum))

##plot it
annotation_flow_plot <- ggplot(lodes, aes(x = x, y = n_cells, alluvium = alluvium, stratum = stratum)) +
  geom_alluvium() +
  geom_stratum(aes(fill = after_stat(stratum))) +
  geom_text(stat = "stratum", aes(label = after_stat(str_remove(stratum, "\\.\\d+$")))) +
  scale_fill_manual(values = lode_colors) +
  labs(x = NULL, y = "Number of cells") +
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#plot annotations per day
annotations_per_day_plot <- seuratObj@meta.data %>%
  select(pbm_day, assigned_morphotype_final)%>%
  group_by(pbm_day, assigned_morphotype_final)%>%
  summarise(ncells = n(), .groups = "drop") %>%
  complete(pbm_day, assigned_morphotype_final, fill = list(ncells = 0)) %>%
  group_by(pbm_day) %>%
  mutate(percent = ncells/sum(ncells)) %>%
  mutate(percent_text = paste0(round(percent*100, 1), "%"))%>%
  ggplot(aes(x = pbm_day, alluvium = assigned_morphotype_final, stratum = assigned_morphotype_final, y = percent, fill = assigned_morphotype_final, label = percent_text))+
  geom_alluvium(alpha = 0.8, decreasing = FALSE)+geom_text_repel(stat = "stratum", decreasing = FALSE, direction = "y", nudge_x = 0.1)+
  labs(fill = "Morphotype", x = "Library", y = "Percentage")+
  scale_fill_manual(values = seuratObj@misc$colors$assigned_morphotype_final)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#format the umap plots
umaps <- lapply(umaps, function(x){
  x <- x +
    scale_x_continuous(breaks = NULL, name = "UMAP 1") +
    scale_y_continuous(breaks = NULL, name = "UMAP 2") +
  return(x)
})


#now build the figure panel
plot_list <- list()
plot_list[["A"]] <- ggalign::align_plots(umaps[[1]], free_lab(umaps[[2]]), ncol = 1) + layout_tags(NULL)
plot_list[["B"]] <- cluster_mapping_plot
plot_list[["C"]] <- annotation_flow_plot
plot_list[["D"]] <- umaps[[3]]
plot_list[["E"]] <- annotations_per_day_plot
right <- align_plots(plot_list[["B"]], plot_list[["C"]], ncol = 1) 
top <- align_plots(plot_list[["A"]], right, nrow = 1) 
bottom <- align_plots(plot_list[["D"]], plot_list[["E"]], nrow = 1)

final_plot <- align_plots(top, bottom, ncol = 1, heights = c(0.7, 0.3)) + layout_tags("A")
#final_plot

#save the figure
dir.create("../Figures")
ggsave("../Figures/Figure1.pdf", final_plot, width = 8.27, height = 7, scale = 2)
