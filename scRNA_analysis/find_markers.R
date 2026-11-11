#this script find marker genes for clusters of cells

#prepare the environment####
rm(list = ls(all.names = TRUE))
gc()

#load libraries####
library(glmGamPoi)
library(Seurat)
library(MAST)
library(tidyverse)
library(conflicted)
source("git_modules/scRNA_tools/scRNA_functions.R")

conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("vars", "dplyr")
conflict_prefer("as.factor", "base")
conflict_prefer("align_plots", "ggalign")

#get color pallets
source("color_pallets.R")
#just for readability:
morphotype_colors <- morphotype_colors
cell_cycle_colors <- cell_cycle_colors

#set paths
seuratObj_path <- "../inputs/scRNA/seurat_objects/seuratObj_all_timepoints_processed.rds"
out_dir <- "../outputs/scRNA/Heatmaps"

##set parameters####
categories_to_group <- c("seurat_clusters", "atlas.stages.ident", "assigned_morphotype_final")

#set outputs paths based on input object name####
base_name <- gsub("_unprocessed|.rds", "", basename(seuratObj_path))
out_seurat_name <- paste0(base_name, "_processed.rds")

##fix issues with memory limit####
library(future)
options(future.globals.maxSize = (32*(1024^3)))

##create output directory if needed####
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#actual analysis####
##read the seurat object####
seuratObj <- readRDS(file = seuratObj_path)

#temporary: set cluster 9 to `gamete-like`
seuratObj@meta.data <- seuratObj@meta.data %>%
  mutate(assigned_morphotype_final = ifelse(seurat_clusters == "9", "gamete-like", assigned_morphotype_final))

#run SCT on it
seuratObj[["RNA"]] <- split(seuratObj[["RNA"]], f = seuratObj$orig.ident)
seuratObj <- SCTransform(seuratObj, return.only.var.genes = FALSE)
seuratObj <- PrepSCTFindMarkers(seuratObj)

#check marker genes####
for(Category in categories_to_group){
  marker_genes <- FindAllMarkers(seuratObj, group.by = Category, logfc.threshold = log2(2), test.use = "MAST", only.pos = T, min.pct = 0.1)
  marker_genes <- marker_genes %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    as.data.frame() %>%
    mutate(transcript_id = gene) %>%
    dplyr::rename(gene_id = gene) %>%
    mutate(gene_id = gsub("LdCL-", "LdCL_", gene_id)) %>%
    mutate(gene_id = gsub("LdBPK-", "LdBPK_", gene_id)) %>%
    mutate(gene_id = gsub("-t42-1", "", gene_id)) %>%
    mutate(trypDB_link = ifelse(grepl("LdCL", gene_id), paste0("https://tritrypdb.org/tritrypdb/app/record/gene/", gene_id), NA)) %>%
    left_join(seuratObj[["RNA"]]@meta.data %>% select(-transcript_id) %>% rownames_to_column("transcript_id"), by = "transcript_id", keep =  FALSE)
  
  #save excel file with marker_genes
  xlsx::write.xlsx(marker_genes, paste0(out_dir, "/", Category, "_markers.xlsx"))
  
  genes_to_plot <- marker_genes %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 20) %>%
    pull("transcript_id") %>%
    unique()
  
  #plot heatmap####
  gene_names <- seuratObj@assays$RNA@meta.data %>%
    rownames_to_column("rowname") %>%
    select(rowname, cellranger_name) %>%
    deframe()
  
  plot <- DoHeatmap_plus(seuratObj, 
                         features = genes_to_plot,
                         feature_labels = gene_names,
                         group.by = c(Category, categories_to_group[categories_to_group != Category]),
                         group.colors = seuratObj@misc[["colors"]]
  )
  
  
  ggsave(plot = plot, filename = paste0(out_dir, "/", Category, "_markers_Heatmap.pdf"), height = 20, width = 30)
}
