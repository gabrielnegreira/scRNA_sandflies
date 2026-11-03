#this script does the general analysis, like dimensionality reduction, DE genes, etc.

#prepare the environment####
rm(list = ls(all.names = TRUE))

##load libraries####
#library(glmGamPoi)
library(Seurat)
library(MAST)
library(tidyverse)
library(clustree)
library(Cairo)
library(ggalluvial)
library(ggrepel)
library(ggalign)
source("git_modules/scRNA_tools/scRNA_functions.R")
source("git_modules/GO_tools/GO_tools.R")

#get paths to external scripts and sources
source("scRNA_analysis/env.R")

##set parameters####
resolutions_to_try <- c(0:10)/10
chosen_resolution <- 0.4
pca_dims <- 1:20

#set paths
#seuratObj_path <- "../inputs/scRNA/seurat_objects/B12_E4_G3_G6_seurat_unprocessed.rds"
seuratObj_path <- "../inputs/scRNA/seurat_objects/B12_E4_seurat_unprocessed.rds"
out_dir <- "../outputs/scRNA/"

#set outputs paths based on input object name
out_seurat_name <- gsub("unprocessed", "processed", basename(seuratObj_path))
out_dir <- paste0(out_dir, gsub("_unprocessed.rds", "", basename(seuratObj_path)))

##fix issues with memory limit####
library(future)
options(future.globals.maxSize = (32*(1024^3)))

##create output directory if needed####
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#actual analysis####
##read the seurat object####
seuratObj <- readRDS(file = seuratObj_path)

##split the seurat object per library####
#first part of the analysis must be done on each library individually, so here we split the seuratobj by library
misc_data <- seuratObj@misc #this is lost after splitting so we save it before.
seuratObj <- SplitObject(seuratObj, split.by = "orig.ident")

##analyse each library separated####

###remove outlier cells and lowly expressed genes####
#####update cells metadata to include expression of mitochondrial and ribosomal genes
mito_features <- rownames(filter(seuratObj[[1]][["RNA"]]@meta.data, is_mitochondrial))
ribo_features <- rownames(filter(seuratObj[[1]][["RNA"]]@meta.data, is_ribosomal))

seuratObj <- lapply(seuratObj, PercentageFeatureSet, features = mito_features , col.name = "percent_mito")
seuratObj <- lapply(seuratObj, PercentageFeatureSet, features = ribo_features , col.name = "percent_ribo")


###flag outlier cells and remove them
seuratObj <- lapply(seuratObj, function(x){
  x@meta.data <- x@meta.data %>%
    rownames_to_column("cell_barcode") %>%
    mutate(is_outlier_mito = percent_mito > quantile(percent_mito, 0.99)) %>%
    mutate(is_outlier_counts = nCount_RNA > quantile(nCount_RNA, 0.99)) %>%
    mutate(is_outlier_features = nFeature_RNA > quantile(nFeature_RNA, 0.99)) %>%
    ungroup() %>%
    column_to_rownames("cell_barcode")
  x@meta.data$to_remove <- rowSums(x@meta.data[,grep("is_outlier_", colnames(x@meta.data))]) >= 1
  x <- subset(x, !to_remove)
  return(x)
})

###Normalize and scale genes####
seuratObj <- seuratObj %>%
  lapply(NormalizeData) %>%
  lapply(ScaleData) %>%
  lapply(SCTransform, return.only.var.genes = FALSE)

##join the datasets####
seuratObj <- merge(seuratObj[[1]], seuratObj[-1])
seuratObj@misc <- misc_data #restore misc data
rownames(seuratObj[["RNA"]]@meta.data) <- rownames(seuratObj[["RNA"]]) #restore rownames in metadata (lost after split/merge)

##find variable features####
seuratObj <- FindVariableFeatures(seuratObj)

##fix per-library SCT####
seuratObj <- PrepSCTFindMarkers(seuratObj)

##run PCA####
seuratObj <- RunPCA(seuratObj, verbose = FALSE)
##choose number of PCAs to use for downstream analysis
options <- c(50,  max(pca_dims) + 20)
plot <- ElbowPlot(seuratObj, ndims = options[which.min(options)])+
  geom_vline(xintercept = max(pca_dims), linetype = "dashed", color = "darkred")

CairoPDF(paste0(out_dir, "/elbow_plot.pdf"), height = 6, width = 8)
show(plot)
dev.off()
show(plot)

#cluster cells####
set.seed(123)
seuratObj <- FindNeighbors(seuratObj, dims = pca_dims, nn.method = "annoy") %>%
  FindClusters(resolution = resolutions_to_try)

###choose resolution
plot <- clustree(seuratObj)
plot <- plot+
  geom_hline(yintercept = length(resolutions_to_try) - which(resolutions_to_try == chosen_resolution), linetype = "dashed")+
  ggtitle(paste0("Chosen resolution: ", chosen_resolution, " (dashed line)"))

CairoPDF(paste0(out_dir, "/clustree_plot.pdf"), height = 10, width = 8)
show(plot)
dev.off()
show(plot)


seuratObj$seurat_clusters <- seuratObj[[paste0("SCT_snn_res.",chosen_resolution)]]
seuratObj <- SetIdent(seuratObj, value = "seurat_clusters")
seuratObj@misc$colors$seurat_clusters <- create_colors(seuratObj$seurat_clusters)

#Run UMAP####
seuratObj <- RunUMAP(seuratObj, dims = pca_dims, min.dist = 0.3, n.neighbors = 30)
plot <- DimPlot(seuratObj, split.by = "subclone", cols = seuratObj@misc$colors$seurat_clusters)
CairoPDF(paste0(out_dir, "/UMAP_clusters.pdf"), height = 4, width = 16)
show(plot)
dev.off()
show(plot)

#Run PHATE####
seuratObj <- RunPHATE(seuratObj, dims = pca_dims)
plot <- DimPlot(seuratObj, split.by = "subclone", cols = seuratObj@misc$colors$seurat_clusters)
CairoPDF(paste0(out_dir, "/UMAP_clusters.pdf"), height = 4, width = 16)
show(plot)
dev.off()
show(plot)

#define cell cycle
source(cell_cycle_script)

#visualize the proportion of clusters in each clone####
plot <- seuratObj@meta.data %>%
  select(subclone, seurat_clusters)%>%
  group_by(subclone, seurat_clusters)%>%
  summarise(ncells = n())%>%
  group_by(subclone)%>%
  mutate(percent = ncells/sum(ncells)) %>%
  mutate(percent_text = paste0(round(percent*100, 1), "%"))%>%
  ggplot(aes(x = subclone, alluvium = seurat_clusters, stratum = seurat_clusters, y = percent, fill = seurat_clusters, label = percent_text))+
  geom_alluvium(alpha = 0.8, decreasing = FALSE)+geom_text_repel(stat = "stratum", decreasing = FALSE, direction = "y", nudge_x = 0.1)+
  #geom_col()+geom_text(position = position_stack(vjust = 0.5))+
  labs(fill = "Cluster", x = "Sample", y = "Percentage")+
  scale_fill_manual(values = seuratObj@misc$colors$seurat_clusters)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
CairoPDF(paste0(out_dir, "/cluster_proportions.pdf"), height = 6, width = 8)
show(plot)
dev.off()
show(plot)

#visualize the proportion of stages in each clone####
plot <- seuratObj@meta.data %>%
  select(subclone, atlas.stages.ident)%>%
  group_by(subclone, atlas.stages.ident)%>%
  summarise(ncells = n())%>%
  group_by(subclone)%>%
  mutate(percent = ncells/sum(ncells)) %>%
  mutate(percent_text = paste0(round(percent*100, 1), "%"))%>%
  ggplot(aes(x = subclone, alluvium = atlas.stages.ident, stratum = atlas.stages.ident, y = percent, fill = atlas.stages.ident, label = percent_text))+
  geom_alluvium(alpha = 0.8, decreasing = FALSE)+geom_text_repel(stat = "stratum", decreasing = FALSE, direction = "y", nudge_x = 0.1)+
  #geom_col()+geom_text(position = position_stack(vjust = 0.5))+
  labs(fill = "Cluster", x = "Sample", y = "Percentage")+
  scale_fill_manual(values = seuratObj@misc$colors$atlas.stages.ident)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
CairoPDF(paste0(out_dir, "/morphotype_proportions.pdf"), height = 6, width = 8)
show(plot)
dev.off()
show(plot)

#visualize the proportion of aneuploid cells in each cluster
plot <- seuratObj@meta.data %>%
  group_by(subclone, seurat_clusters, aneu_group) %>%
  summarise(ncells = n()) %>%
  group_by(subclone, seurat_clusters) %>%
  mutate(prop = ncells/sum(ncells)) %>%
  ggplot(aes(x = seurat_clusters, y = ncells, fill = aneu_group))+
  geom_col()+
  facet_wrap(vars(subclone), scales = "free")+
  scale_fill_manual(values = seuratObj@misc$colors$aneu_group)
CairoPDF(paste0(out_dir, "/aneuploid_proportions.pdf"), height = 6, width = 8)
show(plot)
dev.off()
show(plot)

#visualize the proportion of cells in each cell cycle phase in each cluster
plot <- seuratObj@meta.data %>%
  group_by(subclone, seurat_clusters, aneu_group, cell_cycle_phase) %>%
  summarise(ncells = n()) %>%
  group_by(subclone, seurat_clusters, aneu_group) %>%
  mutate(prop = ncells/sum(ncells)) %>%
  ggplot(aes(x = seurat_clusters, y = ncells, fill = cell_cycle_phase))+
  geom_col() +
  facet_wrap(vars(subclone), scales = "free")+
  scale_fill_manual(values = seuratObj@misc$colors$cell_cycle_phase)
CairoPDF(paste0(out_dir, "/cell_cycle_proportions.pdf"), height = 6, width = 8)
show(plot)
dev.off()
show(plot)

#enrichment analysis samples in each cluster####
##create a key/value pair setting the groupings for enrichment analysis
enrich_groups <- c(seurat_clusters = "subclone",
                   seurat_clusters = "aneu_group",
                   subclone = "seurat_clusters",
                   aneu_group = "seurat_clusters")

plot_list <- list()
enrichment_results_list <- list()
for(i in seq_along(enrich_groups)){
  group_a <- names(enrich_groups[i])
  group_b <- unname(enrich_groups[i])
  comparison <- paste(group_a, "Vs", group_b)
  
  ref <- seuratObj@meta.data[[group_a]]
  enrichment_results <- list()
  
  for(group_b_level in unique(seuratObj@meta.data[[group_b]])){
    set <- seuratObj@meta.data %>%
      filter(!!sym(group_b) == group_b_level) %>%
      pull(group_a)
    enrichment <- GO_enrich(set_terms = set, ref_terms = ref, test_type = "both") %>%
      dplyr::rename(!!sym(group_a) := term) %>%
      dplyr::mutate(!!sym(group_b) := group_b_level)
    
    enrichment_results[[group_b_level]] <- enrichment
  }
  
  enrichment_results <- bind_rows(enrichment_results)
  enrichment_results_list[[comparison]] <- enrichment_results
  
  #make the plot
  plot_list[[comparison]] <- enrichment_results %>%
    ggplot(aes(x = log2enrichment, y = .data[[group_a]], color = pvalue < 0.05))+
    geom_segment(aes(x = 0, xend = log2enrichment))+
    geom_point(aes(size = -log10(pvalue)))+
    scale_color_manual(breaks = c("TRUE", "FALSE"), values = c(`TRUE` = "black", `FALSE` = "grey"), limits = c("TRUE", "FALSE"), drop = FALSE)+
    #guides(color = guide_legend(override.aes = list(linetype = 1, linewidth = 0.6, shape = 16)))+ #this forces a line and a dot to be shown in the legend even if there are no `FALSE` values in the plot
    facet_wrap(vars(.data[[group_b]]), ncol = 1)+
    labs(x = "Log2 Enrichment", y = group_a)+
    ggtitle(comparison)
}

plot <- ggalign::align_plots(!!!plot_list, nrow = 2, heights = c(0.3, 0.7), guides = "r") + layout_tags("A")
ggsave(paste0(out_dir, "/subclone_enrichment_by_cluster.pdf"), plot = plot, height = 20, width = 10)

#now do enrichment analysis per clone
enrichment_results <- list()
for(Subclone in unique(seuratObj@meta.data$subclone)){
  ref <- seuratObj@meta.data %>%
    filter(subclone == Subclone) %>%
    pull(aneu_group)
  
  enrichment_results[[Subclone]] <- list()
  for(Cluster in unique(seuratObj@meta.data$seurat_clusters)){
    set <- seuratObj@meta.data %>%
      filter(subclone == Subclone & seurat_clusters == Cluster) %>%
      pull(aneu_group)
    enrichment <- GO_enrich(set, ref, test_type = "both")
    enrichment$cluster <- Cluster
    enrichment$subclone <- Subclone
    enrichment_results[[Subclone]][[Cluster]] <- enrichment
  }
}
enrichment_results <- lapply(enrichment_results, bind_rows)
enrichment_results <- bind_rows(enrichment_results)

#plot it 
x_limits <- ceiling(max(abs(enrichment_results$log2enrichment)))
x_limits <- c(-x_limits, x_limits)
plot <- enrichment_results %>%
  ggplot(aes(x = log2enrichment, y = cluster, color = pvalue < 0.5))+
  geom_point(aes(size = -log10(pvalue)))+
  geom_segment(aes(x = 0, xend = log2enrichment))+
  scale_x_continuous(limits = x_limits)+
  scale_color_manual(breaks = c("TRUE", "FALSE"), values = c(`TRUE` = "black", `FALSE` = "grey"), limits = c("TRUE", "FALSE"), drop = FALSE)+
  facet_grid(cols = vars(subclone), rows = vars(term))

ggsave(filename = paste0(out_dir, "/aneuploidy_enrichment_per_cluster.pdf"), plot = plot, height = 10, width = 16)

#check marker genes####
marker_genes <- FindAllMarkers(seuratObj, logfc.threshold = log2(2), test.use = "MAST", only.pos = T, min.pct = 0.1)
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
xlsx::write.xlsx(marker_genes, paste0(out_dir, "/cluster_markers.xlsx"))

#store marker genes into seurat Object
seuratObj@misc$marker_genes <- marker_genes

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
                       group.by = c("seurat_clusters", "aneu_group", "subclone"),
                       group.colors = seuratObj@misc$colors
                       )


ggsave(plot = plot, filename = paste0(out_dir, "/Heatmap_marker_genes.pdf"), height = 20, width = 30)

#perform GO analysis####
GO_terms <- read_delim(ref_GO_annotations) %>%
  mutate(gene_id = gsub("_", "-", gene_id)) %>%
  mutate(gene_id = paste0(gene_id, "-t42-1"))

GO_meta <- get_GO_data(GO_terms$simplified_term) %>%
  select(input_term, name, aspect)

#Append the metadata
GO_terms <- GO_terms %>%
  mutate(name = GO_meta[simplified_term,]$name,
         aspect = GO_meta[simplified_term,]$aspect)

#now check enriched terms for each cluster
GO_results <- c()
for(go_type in unique(GO_terms$aspect)){
  print(go_type)
  if(is.na(go_type)){
    next
  }
  ref_terms <- filter(GO_terms, aspect == go_type)
  for(cluster_num in unique(marker_genes$cluster)){
    gene_ids <- marker_genes %>%
      filter(cluster == cluster_num & transcript_id %in% ref_terms$gene_id) %>%
      slice_max(order_by = avg_log2FC, n = 50) %>%
      pull(transcript_id)
    if(length(gene_ids) == 0){
      next
    }else{
      set_terms <- ref_terms %>% 
        filter(gene_id %in% gene_ids) %>%
        pull(simplified_term)
      GO_enr_res <- GO_enrich(set_terms, ref_terms$simplified_term, test_type = "enrichment")
      GO_enr_res$cluster <- cluster_num
      GO_enr_res$aspect <- go_type
      GO_results <- rbind(GO_results, GO_enr_res)
    }
  }
}
GO_results$GO_name <- GO_terms$name[match(GO_results$term, GO_terms$simplified_term)]
#quick plot
plot_list <- list()
scale_labels <- format(10^-c(0:4), scientific = FALSE)
scale_labels[5] <- paste0("< ", scale_labels[5])
for(type in unique(GO_results$aspect)){
  plot_list[[type]] <- GO_results %>%
    filter(log2enrichment > 0) %>%
    filter(aspect == type)%>%
    mutate(aspect = gsub("_", " ", aspect)) %>%
    group_by(aspect) %>%
    mutate(GO_name = paste0(GO_name, ": (", term, ")")) %>%
    mutate(cluster = ifelse(!is.na(suppressWarnings(as.numeric(cluster))), as.numeric(cluster), cluster)) %>%
    mutate(GO_name = fct_reorder2(GO_name, dplyr::desc(pvalue), cluster)) %>%
    mutate(pvalue = -log10(pvalue)) %>%
    mutate(pvalue = ifelse(pvalue > 4, 4, pvalue)) %>%
    ggplot(aes(x = as.factor(cluster), y = GO_name, fill = pvalue))+
    geom_tile()+
    labs(fill = "p-value", x = "cell cluster")+
    guides(fill = guide_colorbar(reverse = TRUE))+
    facet_grid(rows = vars(aspect), scales = "free", space = "free")+
    scale_fill_gradientn(colors = c("white", "purple"), limits = c(0, 4), breaks = c(0:4), labels = scale_labels)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(panel.background = element_blank(), axis.title.y = element_blank())
}
panel_props <- unlist(lapply(plot_list, function(x)length(unique(x$data$term))))
panel_props <- sqrt(panel_props)
panel_props <- panel_props/sum(panel_props)
CairoPDF(paste0(out_dir, "/marker_genes_GO.pdf"), width = 9, height = 18)
show(plot_grid(plotlist = plot_list, ncol = 1, align = "vh", rel_heights = panel_props))
dev.off()

#save the seurat object to disk
saveRDS(seuratObj, file = paste0(dirname(seuratObj_path), "/", out_seurat_name))