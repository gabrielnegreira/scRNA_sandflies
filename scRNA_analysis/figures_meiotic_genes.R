#clean environment
source("clean_environment.R")

#library####
library(tidyverse)
library(Seurat)
library(ggalign)
source("custom_functions.R")

#parameters####
seuratObj <- readRDS("project_specific/BPK081_subclones/outputs/all_subclones_20_dims_annoy_7_clusters/seuratObj.rds")
out_dir <- "project_specific/BPK081_subclones/"
meiotic_genes_LTR <- c(
  #HOP1    = "LTRL590_360017200",
  #SPO11   = "LTRL590_160011700",
  #DMC1    = "LTRL590_350055700",
  #HOP2    = "LTRL590_270031700",
  #MND1    = "LTRL590_240015900",
  `HAP2-1` = "LTRL590_350009500",
  `HAP2-2` = "LTRL590_360048900"
  )

LTR_to_LdCL <- read_delim("common_inputs/LdCL_ref/LTR_to_LdCL.txt") %>%
  dplyr::rename(LdCL = source_id, 
                LTR = `Input Ortholog(s)`, 
                ortholog_group = `Ortholog Group`) %>%
  dplyr::select(LdCL, LTR, ortholog_group) %>%
  mutate(LdCL = gsub("_", "-", LdCL)) %>%
  separate_longer_delim(cols = LTR, delim = ",") %>%
  select(LTR, LdCL) %>%
  deframe()

meiotic_genes <- LTR_to_LdCL[meiotic_genes_LTR]
names(meiotic_genes) <- names(meiotic_genes_LTR)


#if not done, join layers
seuratObj <- JoinLayers(seuratObj, assay = "RNA")

#get the data needed for plotting
umap <- Embeddings(seuratObj, reduction = "umap") %>%
  as_tibble(rownames = "cell")

cell_meta <- seuratObj@meta.data %>%
  select(subclone, aneu_group, seurat_clusters) %>%
  as_tibble(rownames = "cell")

expression <- GetAssayData(seuratObj, assay = "RNA", layer = "data")[meiotic_genes,] %>%
  t() %>%
  as_tibble(rownames = "cell") %>%
  pivot_longer(cols = -cell, names_to = "gene", values_to = "expression")

#combine_them
to_plot <- expression %>%
  left_join(umap, by = "cell") %>%
  left_join(cell_meta, by = "cell")

#add gene names
gene_names <- names(meiotic_genes)
names(gene_names) <- meiotic_genes
to_plot$gene_name <- gene_names[to_plot$gene]

#calculate cluster centers in the umap space to add the labels
cluster_labels <- to_plot %>%
  group_by(subclone, seurat_clusters) %>%
  summarise(umap_1 = density_peak(umap_1), umap_2 = density_peak(umap_2))

#plot 
plot_list <- list()
plot_list[["umap_clusters"]] <- to_plot %>%
  ggplot(aes(x = umap_1, y = umap_2, color = seurat_clusters))+
  geom_point(size = 0.2)+
  geom_text(data = cluster_labels, aes(label = seurat_clusters), color = "black", size = 5)+
  scale_color_manual(values = seuratObj@misc$colors$seurat_clusters)+
  facet_grid(cols = vars(subclone))+
  theme_void()+
  labs(color = "Cluster")+
  ggtitle("Seurat Clusters")+
  #guides(color = guide_legend(override.aes = list(size = 2)))+
  guides(color = "none")+
  theme(strip.text = element_text(size = 16))

plot_list[["umap_aneuploidy"]] <- to_plot %>%
  ggplot(aes(x = umap_1, y = umap_2, color = aneu_group))+
  geom_point(size = 0.2)+
  scale_color_manual(values = seuratObj@misc$colors$aneu_group)+
  facet_grid(cols = vars(subclone))+
  theme_void()+
  labs(color = "aneu_group")+
  ggtitle("Aneuploidy groups")+
  guides(color = guide_legend(override.aes = list(size = 3)))+
  theme(strip.text = element_blank())


plot_list[["umap_expression"]] <- to_plot %>%
  group_by(gene) %>%
  mutate(expression = ifelse(expression > quantile(expression, 0.99), quantile(expression, 0.99), expression)) %>% #cap the scale at the 99% percentil
  mutate(expression = expression/max(expression)) %>% #normalize expression per gene
  ungroup() %>%
  arrange(subclone, gene, expression) %>%
  mutate(expression = ifelse(expression == 0, NA, expression)) %>%
  ggplot(aes(x = umap_1, y = umap_2, color = expression))+
  geom_point(size = 0.2)+
  theme_void()+
  facet_grid(cols = vars(subclone), rows = vars(gene_name), scales = "free", switch = "y")+
  labs(color = "Relative\nExpression")+
  scale_color_gradientn(colors = c("#fff9cc", "#fde954", "#ea5436", "#445ea1", "#20294f"), na.value = "lightgrey")+
  ggtitle("Expression")+
  theme(strip.text.x.top = element_blank(),
        strip.text.y.left = element_text(angle = 90, size = 16))

#now plot the violin 
plot_list[["violin_subclones"]] <- to_plot %>%
  group_by(gene) %>%
  mutate(expression = expression + (rnorm(n()) / 1e5)) %>% #adds a small noise to the data, same as seurat does
  ggplot(aes(x = subclone, y = expression, group = paste0(subclone, aneu_group), fill = aneu_group))+
  geom_violin(scale = "width")+
  geom_point(size = 0.01, alpha = 0.5, position = position_jitterdodge(dodge.width = 0.9), show.legend = FALSE)+
  facet_grid(rows = vars(gene_name))+
  scale_fill_manual(values = seuratObj@misc$colors$aneu_group)+
  theme_cowplot()+
  labs(fill = NULL)+
  ggtitle("Expression by aneu_group")+
  theme(strip.text = element_text(hjust = 0, size = 16))

plot_list[["violin_clusters"]] <- to_plot %>%
  group_by(gene) %>%
  mutate(expression = expression + (rnorm(n()) / 1e5)) %>% #adds a small noise to the data, same as seurat does
  ggplot(aes(x = subclone, y = expression, group = paste(seurat_clusters, subclone), fill = seurat_clusters))+
  geom_violin(scale = "width")+
  geom_point(size = 0.01, alpha = 0.5, position = position_jitterdodge(dodge.width = 0.9), show.legend = FALSE)+
  facet_grid(rows = vars(gene_name))+
  scale_fill_manual(values = seuratObj@misc$colors$seurat_clusters)+
  theme_cowplot()+
  labs(fill = "Cluster")+
  ggtitle("Expression by cluster")+
  theme(strip.text = element_text(hjust = 0, size = 16))




#adjust text size
plot_list <- lapply(plot_list, function(x){
  x + theme(text = element_text(size = 20))
})

#make final assembly
nrows <- 2 + length(meiotic_genes)
panel_heigths <- c(1/nrows, 1/nrows, length(meiotic_genes)/nrows)
left <- ggalign::align_plots(!!!plot_list[1:3], ncol = 1, heights = panel_heigths) + layout_tags(NULL)
right <- ggalign::align_plots(plot_list$violin_subclones, NULL, plot_list$violin_clusters, ncol = 1, heights = c(0.45, 0.1, 0.45))
#right <- ggalign::align_plots(NULL, right, NULL, ncol = 1, heights = c(0.25, 0.5, 0.25))
right <- ggalign::align_plots(NULL, right, nrow = 1, widths = c(0.1, 0.9))
final <- ggalign::align_plots(left, right, nrow = 1, widths = c(0.6, 0.4)) + 
  layout_tags("A") +
  layout_theme(plot.tag = element_text(size = 20))

#save the final figure panel####
plot_scale <- 4
ggsave(paste0(out_dir, "/meiotic_genes_umap.pdf"), plot = final, width = 8.27 * plot_scale, height = 4 * plot_scale)

#open it
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", args = paste0(out_dir, "/meiotic_genes_umap.pdf"), wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
  system2("xdg-open", args = paste0(out_dir, "/meiotic_genes_umap.pdf"), wait = FALSE)
}

