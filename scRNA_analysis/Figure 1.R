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
seuratObj <- readRDS("../inputs/scRNA/seurat_objects/seuratObj_all_timepoints_processed.rds")

#temporary: test different cluster resolutions
#seuratObj$seurat_clusters <- seuratObj$RNA_snn_res.0.5
#seuratObj@misc$colors$seurat_clusters <- create_colors(levels(seuratObj$seurat_clusters))

#temporary: set a custom color for pbm_day
seuratObj@misc$colors$pbm_day <- pbm_day_colors

#plot cluster umap####
umaps <- list()
umaps[[1]] <- DimPlot(seuratObj, group.by = "pbm_day", cols = seuratObj@misc$colors$pbm_day, label = TRUE) + ggtitle("Clusters")

#plot umap per day
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
  guides(name = "Cluster", color = guide_legend(override.aes = list(size = 2.5), keyheight = grid::unit(0.25, "cm"), keywidth = grid::unit(0.25, "cm")))+
  theme_void()+
  theme(legend.spacing.y = grid::unit(0, "cm"), legend.key.height = grid::unit(0.25, "cm"))

umaps[[2]] <- DimPlot(seuratObj, group.by = "atlas.stages.ident", cols = seuratObj@misc$colors$atlas.stages.ident) + ggtitle("Raw Mapping to Atlas")

#plot the morphotype_to_cluster mapping
cluster_mapping_plot <- seuratObj@meta.data %>%
  group_by(seurat_clusters, assigned_morphotype_final) %>%
  summarise(n_cells = n()) %>%
  ggplot(aes(x = seurat_clusters, y = n_cells, fill = assigned_morphotype_final))+
  geom_col()+
  scale_fill_manual(values = seuratObj@misc$colors$atlas.stages.ident, name = "Morphotype")

#plot umap with final assigned morphotypes
umaps[[3]] <- DimPlot(seuratObj, group.by = "assigned_morphotype_final", cols = seuratObj@misc$colors$assigned_morphotype_final)

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
  mutate(percent_text = ifelse(percent < 0.1, NA, percent_text)) %>%
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
    ggtitle(NULL)
  return(x)
})

#reduce the legend in other plots
plot_list <- list(cluster_mapping_plot, umaps[[3]], annotations_per_day_plot)
plot_list <- lapply(plot_list, function(x){
  x <- x + theme(
    legend.key.size = unit(0.5, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.text = element_text(size = 6)
  )
  return(x)
})

cluster_mapping_plot <- plot_list[[1]]
umaps[[3]] <- plot_list[[2]]
annotations_per_day_plot <- plot_list[[3]]

#save supplementary figures (for later)
panel_list <- list()
panel_list[["S1"]] <- umaps[[2]]
panel_list[["S2"]] <- annotation_flow_plot

#now build the main figure panel

panel_list[["A"]] <- umaps[[1]]
panel_list[["B"]] <- umap_per_day
panel_list[["C"]] <- cluster_mapping_plot
panel_list[["D"]] <- umaps[[3]]
panel_list[["E"]] <- annotations_per_day_plot
panel_list <- lapply(panel_list, free_lab)
top_right_bottom <- align_plots(free_lab(panel_list[["C"]]), panel_list[["D"]], panel_list[["E"]], nrow = 1) 
top_right <- align_plots(panel_list[["B"]], top_right_bottom, ncol = 1, heights = c(0.4, 0.6))
top <- align_plots(panel_list[["A"]], top_right, nrow = 1, widths = c(0.3, 0.7))
bottom <- align_plots(panel_list[["D"]], panel_list[["E"]], nrow = 1)

#final_plot <- align_plots(top, bottom, ncol = 1, heights = c(0.7, 0.3)) + layout_tags("A")
final_plot <- top

#save the figure
dir.create("../Figures")
ggsave("../Figures/Figure1.pdf", final_plot, width = 8.27, height = 3, scale = 2)
