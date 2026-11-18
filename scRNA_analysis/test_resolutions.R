#this is an auxiliary script meant to test different resolutions
#this script find marker genes for morphotypes and seurat clusters

#prepare the environment####
rm(list = ls(all.names = TRUE))
gc()

#load libraries####
library(Seurat)
library(glmGamPoi)
library(presto)
library(MAST)
library(tidyverse)
library(clustree)
library(conflicted)
source("git_modules/scRNA_tools/scRNA_functions.R")

#solve package conflicts####
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("vars", "dplyr")
conflict_prefer("as.factor", "base")
conflict_prefer("align_plots", "ggalign")

#set parameters####
resolutions_to_try <- seq(0, 1, by = 0.05)
out_dir <- "../outputs/scRNA/testing_different_resolutions/"
seuratObj <- readRDS("../inputs/scRNA/seurat_objects/sandfly_data/sandfly_processed.rds")


#do different clusters####
seuratObj <- FindClusters(seuratObj, resolution = resolutions_to_try)

#save clustree####
clustree_plot <- clustree(seuratObj)
ggsave(paste0(out_dir, "/clustree_plot.pdf"), clustree_plot, width = 10, height = 10)

#get col names of the resolutions in the metadata####
resolutions <- grep(paste0(DefaultAssay(seuratObj), "_snn_res"), colnames(seuratObj@meta.data), value = TRUE)

#plot UMAPs
plot_list <- list()
cluster_colors <- list()
for(Resolution in resolutions){
  cluster_colors[[Resolution]] <- create_colors(unique(seuratObj@meta.data[[Resolution]]))
  plot_list[[Resolution]] <- DimPlot(seuratObj, group.by = Resolution, cols = cluster_colors[[Resolution]], label = TRUE)
}
final_plot <- align_plots(!!!plot_list, ncol =1)
dir.create(out_dir, showWarnings = FALSE)
ggsave(paste0(out_dir, "/different_resolutions_umap.png"), final_plot, width = 8, height = 6 * length(plot_list), dpi = 300, limitsize = FALSE)

#now plot heatmaps####
#check marker genes####
for(Resolution in resolutions[2:length(resolutions)]){
  marker_genes <- FindAllMarkers(
    seuratObj, 
    group.by = Resolution, 
    logfc.threshold = log2(1.5), 
    only.pos = T, 
    min.pct = 0.1
  )
  if(nrow(marker_genes) == 0){
    warning(paste0("No markers found at resolution: ", Resolution))
    next
  }
  marker_genes <- marker_genes %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    as.data.frame() %>%
    mutate(seurat_id = gene) %>%
    dplyr::rename(gene_id = gene) %>%
    mutate(gene_id = gsub("LdCL-", "LdCL_", gene_id)) %>%
    mutate(gene_id = gsub("LdBPK-", "LdBPK_", gene_id)) %>%
    mutate(gene_id = gsub("-t42-1", "", gene_id)) %>%
    mutate(trypDB_link = ifelse(grepl("LdCL", gene_id), paste0("https://tritrypdb.org/tritrypdb/app/record/gene/", gene_id), NA)) %>%
    left_join(seuratObj[["RNA"]]@meta.data, by = "seurat_id", keep =  FALSE)
  
  #save excel file with marker_genes
  xlsx::write.xlsx(marker_genes, paste0(out_dir, "/", Resolution, "_markers.xlsx"))
  
  genes_to_plot <- marker_genes %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 20) %>%
    pull("seurat_id") %>%
    unique()
  
  #plot heatmap####
  gene_names <- seuratObj@assays$RNA@meta.data %>%
    rownames_to_column("rowname") %>%
    select(rowname, annotated_name) %>%
    deframe()
  
  plot <- DoHeatmap_plus(seuratObj, 
                         features = genes_to_plot,
                         feature_labels = gene_names,
                         group.by = Resolution,
                         group.colors = cluster_colors[[Resolution]]
  )
  
  
  ggsave(plot = plot, filename = paste0(out_dir, "/", Resolution, "_markers_Heatmap.pdf"), height = 20, width = 30)
}




