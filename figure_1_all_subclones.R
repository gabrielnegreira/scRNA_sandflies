#library####
library(tidyverse)
library(Seurat)
library(ggalign)
source("GO_tools.R")
source("custom_functions.R")

#parameters####
seuratObj <- readRDS("project_specific/BPK081_subclones/outputs/all_subclones_20_dims_annoy_7_clusters/seuratObj.rds")
#seuratObj <- readRDS("project_specific/BPK081_subclones/outputs/B12_and_E4_only/seuratObj.rds") 
out_dir <- "project_specific/BPK081_subclones/"
GO_terms <- read_delim("common_inputs/LdCL_ref/LdCL_GO_terms_updated.tsv") %>%
  mutate(gene_id = gsub("_", "-", gene_id)) %>%
  mutate(gene_id = paste0(gene_id, "-t42-1"))

#if not done, join layers
seuratObj <- JoinLayers(seuratObj, assay = "RNA")

#plot figures####
##umap figure####
figures <- list()
plot_index <- 1
up <- DimPlot(seuratObj, group.by = "seurat_clusters", split.by = "subclone", cols = seuratObj@misc$colors$seurat_clusters, label = TRUE, pt.size = 0.0001) & 
  theme_void() & 
  ggtitle(NULL) &
  guides(color = "none") &
  theme(legend.title = element_text())

down <- DimPlot(seuratObj, group.by = "aneu_group", split.by = "subclone", cols = seuratObj@misc$colors$aneu_group, pt.size = 0.0001) & 
  theme_void() & 
  ggtitle(NULL) &
  theme(strip.text = element_blank())
figures[[plot_index]] <- ggalign::align_plots(up, down, ncol = 1) + layout_tags(NULL)

##cluster proportions####
up <- seuratObj@meta.data %>%
  select(subclone, seurat_clusters)%>%
  group_by(subclone, seurat_clusters)%>%
  summarise(ncells = n())%>%
  group_by(subclone)%>%
  mutate(percent = ncells/sum(ncells)) %>%
  mutate(percent_text = paste0(round(percent*100, 1), "%"))%>%
  ggplot(aes(x = subclone, y = percent, fill = seurat_clusters, label = percent_text))+
  geom_col()+geom_text(position = position_stack(vjust = 0.5), size = 3)+
  labs(fill = "Cluster", x = "Sample", y = "Percentage")+
  scale_fill_manual(values = seuratObj@misc$colors$seurat_clusters)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##cluster enrichment####
enrichment_results <- list()
ref <- seuratObj@meta.data$aneu_group
for(Cluster in unique(seuratObj@meta.data$seurat_clusters)){
  set <- seuratObj@meta.data %>%
    filter(seurat_clusters == Cluster) %>%
    pull(aneu_group)
  enrichment <- GO_enrich(set, ref, test_type = "both")
  enrichment <- enrichment %>%
    dplyr::rename(aneu_group = term) %>%
    mutate(cluster = Cluster)
  enrichment_results[[Cluster]] <- enrichment
}

enrichment_results <- bind_rows(enrichment_results)
x_limits <- ceiling(max(abs(enrichment_results$log2enrichment)))
x_limits <- c(-x_limits, x_limits)
down <- enrichment_results %>%
  filter(pvalue < 0.05) %>%
  ggplot(aes(x = log2enrichment, y = fct_rev(cluster), color = aneu_group))+
  geom_vline(xintercept = 0, color = "lightgrey", linetype = "dashed")+
  geom_segment(aes(x = 0, xend = log2enrichment))+
  geom_point(aes(size = -log10(pvalue)))+
  scale_x_continuous(limits = x_limits)+
  scale_color_manual(values = seuratObj@misc$colors$aneu_group)+
  labs(x = "Log2 enrichment", y = "Cluster", color = NULL)+
  theme_bw()+
  theme(panel.grid = element_blank())

plot_index <- plot_index + 1
figures[[plot_index]] <- ggalign::align_plots(up, down, ncol = 1)

##heatmap of marker genes####
genes_to_plot <- seuratObj@misc$marker_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  pull("transcript_id") %>%
  unique()

gene_names <- seuratObj@assays$RNA@meta.data %>%
  rownames_to_column("rowname") %>%
  select(rowname, cellranger_name) %>%
  deframe()

plot_index <- plot_index + 1
figures[[plot_index]] <- DoHeatmap_plus(seuratObj, 
                       features = genes_to_plot,
                       feature_labels = gene_names,
                       assay = "SCT",
                       feature_label_size = 5,
                       group.by = c("seurat_clusters", "aneu_group", "subclone"),
                       group.colors = seuratObj@misc$colors)

##plot autophagy genes####
autophagy_genes <- names(gene_names[grep("atg|autophagy", tolower(gene_names))])
seuratObj <- AddModuleScore(seuratObj, list("autophagy_score" = autophagy_genes), name = "autophagy_score")

genes_to_plot <- seuratObj@misc$marker_genes %>%
  filter(cluster == "7") %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 5) %>%
  pull("transcript_id")

to_plot <- seuratObj[["RNA"]]["data"][genes_to_plot, , drop = FALSE] %>%
  as_tibble(rownames = "gene") %>%
  mutate(gene_name = gene_names[gene]) %>%
  pivot_longer(cols = -c(gene, gene_name), names_to = "cell", values_to = "expression") %>%
  left_join(rownames_to_column(seuratObj@meta.data, "cell"), by = "cell")

#plot_index <- plot_index + 1
#figures[[plot_index]] <- to_plot %>%
to_plot %>%
  #filter(seurat_clusters != "7") %>%
  mutate(expression = expression + (rnorm(n()) / 1e5)) %>% #adds a small noise to the data, same as seurat does
  ggplot(aes(x = seurat_clusters, y = expression, fill = aneu_group))+
  geom_violin(scale = "width", adjust = 1)+
  #geom_violin(scale = "count", fill = "white", aes(group = aneu_group))+
  #geom_point(aes(group = aneu_group), size = 0.001, position = position_jitterdodge(dodge.width = 0.9))+
  #geom_quasirandom(size = 1, dodge.width = 0.9, shape = 21, stroke = 0.01)+
  #geom_quasirandom(aes(color = aneu_group), size = 0.6, dodge.width = 0.9, shape = 1, stroke = 0.00001)+
  facet_grid(cols = vars(subclone), rows = vars(gene_name))+
  scale_fill_manual(values = seuratObj@misc$colors$aneu_group)+
  scale_color_manual(values = seuratObj@misc$colors$aneu_group)+
  guides(fill = "none")+
  theme_cowplot()+
  theme(strip.text = element_text(hjust = 0, size = 8))





##GO analysis####
GO_meta <- get_GO_data(GO_terms$simplified_term) %>%
  select(input_term, name, aspect)

#Append the metadata
GO_terms <- GO_terms %>%
  mutate(name = GO_meta[simplified_term,]$name,
         aspect = GO_meta[simplified_term,]$aspect)

#get marker genes
marker_genes <- seuratObj@misc$marker_genes

#now check enriched terms for each cluster
##first focus only on biological processes
ref_terms <- filter(GO_terms, aspect == "biological_process")
GO_results <- c()
for(cluster_num in unique(marker_genes$cluster)){
  #get the id of the marker genes of the cluster
  gene_ids <- marker_genes %>%
    filter(cluster == cluster_num & transcript_id %in% ref_terms$gene_id) %>%
    slice_max(order_by = avg_log2FC, n = 50) %>%
    pull(transcript_id)
  #if no marker gene, go to next cluster
  if(length(gene_ids) == 0){
    next
  }else{
    #get the corresponding GO terms of the markers
    set_terms <- ref_terms %>% 
      filter(gene_id %in% gene_ids) %>%
      pull(simplified_term)
    #perform GO enrichment analysis
    GO_enr_res <- GO_enrich(set_terms, ref_terms$simplified_term, test_type = "enrichment")
    GO_enr_res$cluster <- cluster_num
    GO_enr_res$aspect <- go_type
    GO_results <- rbind(GO_results, GO_enr_res)
  }
}
#add the name of the GO terms (description of what each term mean)
GO_results$GO_name <- GO_terms$name[match(GO_results$term, GO_terms$simplified_term)]

#prepare the labels for pvalues
scale_labels <- format(10^-c(0:4), scientific = FALSE)
scale_labels[5] <- paste0("< ", scale_labels[5])

#plot it
plot_index <- plot_index + 1
 figures[[plot_index]] <- GO_results %>%
    filter(log2enrichment > 0) %>%
    mutate(aspect = gsub("_", " ", aspect)) %>%
    mutate(GO_name = paste0(GO_name, ": (", term, ")")) %>%
    mutate(cluster = ifelse(!is.na(suppressWarnings(as.numeric(cluster))), as.numeric(cluster), cluster)) %>%
    mutate(GO_name = fct_reorder2(GO_name, dplyr::desc(pvalue), cluster)) %>%
    mutate(pvalue = -log10(pvalue)) %>%
    mutate(pvalue = ifelse(pvalue > 4, 4, pvalue)) %>%
    ggplot(aes(x = as.factor(cluster), y = GO_name, fill = pvalue))+
    geom_tile()+
    labs(fill = "p-value", x = "cell cluster")+
    guides(fill = guide_colorbar(reverse = TRUE))+
    scale_fill_gradientn(colors = c("white", "purple"), limits = c(0, 4), breaks = c(0:4), labels = scale_labels)+
    theme(panel.background = element_blank(), axis.title.y = element_blank(), axis.text.y = element_text(size = 6))

##cell_cycle####
plot_index <- plot_index + 1
figures[[plot_index]] <- seuratObj@meta.data %>%
  group_by(seurat_clusters, cell_cycle_phase) %>%
  summarise(ncells = n()) %>%
  group_by(seurat_clusters) %>%
  mutate(prop = ncells/sum(ncells)) %>%
  ggplot(aes(x = seurat_clusters, y = ncells, fill = cell_cycle_phase))+
  geom_col(position = position_fill()) +
  #facet_wrap(vars(subclone), scales = "free")+
  scale_fill_manual(values = seuratObj@misc$colors$cell_cycle_phase)


##make final figure panel####
top <- ggalign::align_plots(free_align(figures[[1]], axes = "l"), figures[[2]], nrow = 1, widths = c(0.8, 0.2))
bottom_right <- ggalign::align_plots(NULL, figures[[4]], ncol = 1)
bottom <- ggalign::align_plots(figures[[3]], bottom_right, nrow = 1, widths = c(0.8, 0.2))

final <- ggalign::align_plots(top, bottom, ncol = 1, heights = c(0.35, 0.65)) +
  layout_tags("A")

#save the final figure panel####
plot_scale <- 1.8
ggsave(paste0(out_dir, "/figure_panel_1_scRNA_all_subclones.pdf"), plot = final, width = 8.27 * plot_scale, height = 7 * plot_scale)

#open it
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", args = paste0(out_dir, "/figure_panel_1_scRNA_all_subclones.pdf"), wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
  system2("xdg-open", args = paste0(out_dir, "/figure_panel_1_scRNA_all_subclones.pdf"), wait = FALSE)
}
