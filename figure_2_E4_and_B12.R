#library####
library(tidyverse)
library(Seurat)
library(ggalign)
library(ggbeeswarm)
source("GO_tools.R")
source("custom_functions.R")

#parameters####
seuratObj <- readRDS("project_specific/BPK081_subclones/outputs/all_subclones_20_dims_annoy_7_clusters/seuratObj.rds") %>%
  subset(subset = subclone %in% c("E4", "B12")) 
out_dir <- "project_specific/BPK081_subclones/"
aneuploidies <- list(B12 = c("Ld01", "Ld05", "Ld10", "Ld13", "Ld20", "Ld23"),
                     E4 = c("Ld01", "Ld13", "Ld20", "Ld23", "Ld26"))

#make sure layers are joint
DefaultAssay(seuratObj) <- "RNA"
seuratObj <- JoinLayers(seuratObj)
DefaultAssay(seuratObj) <- "SCT"

#redo enrichment analysis but now per subclone####
figures <- list()
##cluster enrichment
enrichment_results <- list()
for(Subclone in unique(seuratObj@meta.data$subclone)){
  enrichment_results[[Subclone]] <- list()
  ref <- seuratObj@meta.data %>%
    filter(subclone == Subclone) %>%
    pull(aneu_group)
  
  for(Cluster in unique(seuratObj@meta.data$seurat_clusters)){
    set <- seuratObj@meta.data %>%
      filter(subclone == Subclone) %>%
      filter(seurat_clusters == Cluster) %>%
      pull(aneu_group)
    
    enrichment <- GO_enrich(set, ref, test_type = "both")
    enrichment <- enrichment %>%
      dplyr::rename(aneu_group = term) %>%
      mutate(cluster = Cluster,
             subclone = Subclone)
    enrichment_results[[Subclone]][[Cluster]] <- enrichment
  }  
}

enrichment_results <- enrichment_results %>%
  lapply(bind_rows) %>%
  bind_rows()

x_limits <- ceiling(max(abs(enrichment_results$log2enrichment)))
x_limits <- c(-x_limits, x_limits)

plot_index <- 1
figures[[plot_index]] <- enrichment_results %>%
  filter(pvalue < 0.05) %>%
  ggplot(aes(x = log2enrichment, y = fct_rev(cluster), color = aneu_group))+
  geom_vline(xintercept = 0, color = "lightgrey", linetype = "dashed")+
  geom_segment(aes(x = 0, xend = log2enrichment))+
  geom_point(aes(size = -log10(pvalue)))+
  scale_x_continuous(limits = x_limits)+
  scale_color_manual(values = seuratObj@misc$colors$aneu_group)+
  labs(x = "Log2 enrichment", y = "Cluster", color = NULL)+
  facet_wrap(vars(subclone), nrow = 1)+
  theme_bw()+
  theme(panel.grid = element_blank())

#Redo marker genes per cluster#### 
## For now this will not be used. Maybe for supplementary. 
marker_genes <- FindAllMarkers(seuratObj, 
                               assay = "SCT",
                               latent.vars = c("aneu_group", "subclone"),
                               group.by = "seurat_clusters",
                               logfc.threshold = log2(2),
                               only.pos = TRUE,
                               test.use = "MAST",
                               min.pct = 0.1)

marker_genes <- marker_genes %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

genes_to_plot <- marker_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  pull(gene) %>%
  unique()

gene_names <- seuratObj[["RNA"]]@meta.data %>%
  rownames_to_column("id") %>%
  select(id, cellranger_name) %>%
  deframe()
  
plot_index <- plot_index + 1
figures[[plot_index]] <- DoHeatmap_plus(seuratObj, 
                                        features = genes_to_plot, 
                                        group.by = c("seurat_clusters", "aneu_group", "subclone"), 
                                        feature_labels = gene_names,
                                        feature_label_size = 6,
                                        group.colors = seuratObj@misc$colors)


#look for markers between aneuploid and baseline### 
## For now this will not be used. Maybe for supplementary. 
marker_genes <- FindAllMarkers(seuratObj, 
                               assay = "SCT",
                               latent.vars = c("seurat_clusters", "subclone"),
                               group.by = "aneu_group",
                               logfc.threshold = log2(2),
                               only.pos = TRUE,
                               test.use = "MAST",
                               min.pct = 0.1)

marker_genes <- marker_genes %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC))

genes_to_plot <- marker_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  pull(gene) %>%
  unique()

plot_index <- plot_index + 1
figures[[plot_index]] <- DoHeatmap_plus(seuratObj, 
                                        features = genes_to_plot, 
                                        group.by = c("seurat_clusters", "aneu_group", "subclone"), 
                                        feature_labels = gene_names,
                                        feature_label_size = 6,
                                        group.colors = seuratObj@misc$colors)


##Find DE genes between aneuploids and baselines in each cluster
markers_list <- list()
for(Cluster in unique(seuratObj$seurat_clusters)){
  sub <- subset(seuratObj, seurat_clusters == Cluster)
  DefaultAssay(sub) <- "RNA"
  sub <- JoinLayers(sub)
  DefaultAssay(sub) <- "SCT"
  sub@assays$SCT@SCTModel.list[c(2:length(sub@assays$SCT@SCTModel.list))] <- NULL #fix bug with multiple SCTModels still remaining after `PrepSCTFindMarkers()`.
  markers <-  FindMarkers(
    sub,
    ident.1 = "aneuploid",
    ident.2 = "baseline",
    group.by = "aneu_group",
    logfc.threshold = log2(2),
    min.pct = 0.1,
    assay = "SCT",
    #latent.vars = "subclone",
    test.use = "MAST"   # can be changed to 'MAST' or 'DESeq2' if desired
  )
  markers <- markers %>%
    as.tibble(rownames = "feature") %>%
    mutate(cluster = Cluster)
  markers_list[[Cluster]] <- markers
}

#suspecious markers are those found in aneuploid chromosomes (so they could be simply a consequence of aneuploidy)
markers <- bind_rows(markers_list) %>%
  filter(p_val_adj < 0.05) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  mutate(suspicious = grepl(gsub("Ld", "LdCL-", paste(unique(as.character(unlist(aneuploidies))), collapse = "|")), feature)) %>%
  mutate(suspicious = suspicious & avg_log2FC > 0 & avg_log2FC < log2(3))
  

genes_to_plot <- markers %>%
  #filter(!suspicious) %>%
  #group_by(feature) %>%
  #mutate(occurences = n()) %>%
  #filter(occurences == 1) %>%
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  #slice_head(n = 20) %>%
  pull("feature") %>%
  unique()

plot_index <- plot_index + 1
figures[[plot_index]] <- DoHeatmap_plus(seuratObj, 
                                        features = genes_to_plot, 
                                        group.by = c("seurat_clusters", "aneu_group", "subclone"), 
                                        feature_labels = gene_names,
                                        feature_label_size = 6,
                                        group.colors = seuratObj@misc$colors)

#violin plots of cluster 2 markers####
interesting_genes <- c(
  "LdCL-270017700-t42-1", #hypothetical protein.1786
  "LdCL-160010800-t42-1", #Core histone-like transcription factor
 # "LdCL-010012900-t42-1",  #EIF4a
  "LdCL-290024100-t42-1", #histone H2A putative 4
  "LdCL-340043700-t42-1" #hypothetical protein.2894
  )

to_plot <- seuratObj[["RNA"]]["scale.data"][interesting_genes, , drop = FALSE] %>%
  as_tibble(rownames = "gene") %>%
  mutate(gene_name = gene_names[gene]) %>%
  pivot_longer(cols = -c(gene, gene_name), names_to = "cell", values_to = "expression") %>%
  left_join(rownames_to_column(seuratObj@meta.data, "cell"), by = "cell")

plot_index <- plot_index + 1
figures[[plot_index]] <- to_plot %>%
  filter(seurat_clusters == "2") %>%
  group_by(gene_name) %>%
  mutate(expression = expression + (rnorm(n()) / 1e5)) %>% #adds a small noise to the data, same as seurat does
  ggplot(aes(x = subclone, y = expression, fill = aneu_group))+
  geom_violin(scale = "width")+
  #geom_violin(scale = "count", fill = "white", aes(group = aneu_group))+
  geom_point(aes(group = aneu_group), size = 0.001, position = position_jitterdodge(dodge.width = 0.9))+
  #geom_quasirandom(size = 1, dodge.width = 0.9, shape = 21, stroke = 0.01)+
  #geom_quasirandom(aes(color = aneu_group), size = 0.6, dodge.width = 0.9, shape = 1, stroke = 0.00001)+
  facet_grid(rows = vars(gene_name))+
  scale_fill_manual(values = seuratObj@misc$colors$aneu_group)+
  scale_color_manual(values = seuratObj@misc$colors$aneu_group)+
  guides(fill = "none")+
  theme_cowplot()+
  theme(strip.text = element_text(hjust = 0, size = 8))





#build final plot
left <- ggalign::align_plots(figures[[1]], figures[[5]], ncol = 1, heights = c(0.2, 0.8))
final <- ggalign::align_plots(left, figures[[4]], nrow = 1, widths = c(0.3, 0.7))  + layout_tags("A")

#save the final figure panel####
plot_scale <- 1.8
ggsave(paste0(out_dir, "/figure_panel_2_scRNA_B12_and_E4.pdf"), plot = final, width = 8.27 * plot_scale, height = 5 * plot_scale)

#open it
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", args = paste0(out_dir, "/figure_panel_2_scRNA_B12_and_E4.pdf"), wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
  system2("xdg-open", args = paste0(out_dir, "/figure_panel_2_scRNA_B12_and_E4.pdf"), wait = FALSE)
}
