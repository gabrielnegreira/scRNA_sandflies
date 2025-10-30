#clean environment
source("clean_environment.R")

#library####
source("custom_functions.R")
library(tidyverse)
library(Seurat)
library(DESeq2)
library(ggalign)
library(ggrepel)
library(ggvenn)

#inputs####
seuratObj <- readRDS("project_specific/BPK081_subclones/outputs/all_subclones_20_dims_annoy_7_clusters/seuratObj.rds") #%>%
  #subset(subset = subclone %in% c("B12", "E4"))

#combine G3 and G6 as a single dataset (temporary)####
#seuratObj@meta.data <- seuratObj@meta.data %>%
#  mutate(subclone = ifelse(subclone %in% c("G3", "G6"), "G3_G6", subclone))


#other inputs####
out_dir <- "project_specific/BPK081_subclones/"
aneu_chromo <- list("B12" = c("01", "05", "10", "13", "20", "23"),
                    "E4" = c("13", "20", "23", "26"),
                    "G3" = c("20", "23", "26", "29"),
                    "G6" = c("01", "09", "13", "16", "20", "22", "23", "24", "26", "35"))


#parameters####
group_by <- c("subclone", "seurat_clusters") #first value determines DE group, additional values work as "replicates"
pvalue_threshold <- 0.05
log2FC_threshold <- log2(1.5)
artificial_replicates <- 0
reference_group <- "G3" #this tells which sample is the reference where the other samples will be compared to. 

#set gene names####
gene_names <- seuratObj@assays$RNA@meta.data %>%
  as_tibble(rownames = "gene") %>%
  select(gene, cellranger_name) %>%
  deframe()

#code####
#add artificial replicates
if(artificial_replicates > 1){
  seuratObj$art_rep <- sample(c(1:artificial_replicates), size = ncol(seuratObj), replace = TRUE)
  group_by <- c(group_by, "art_rep")
}

#mark genes which might come from aneuploid chromosomes
to_append <- as.data.frame(matrix(ncol = length(aneu_chromo), nrow = nrow(seuratObj@assays$RNA)))
colnames(to_append) <- names(aneu_chromo)
values <- c("baseline", "aneuploid")
for(Subclone in colnames(to_append)){
  to_append[,Subclone] <- values[((as.integer(seuratObj@assays$RNA@meta.data$chromosome) %in% as.integer(aneu_chromo[[Subclone]]))) + 1]
}

seuratObj@assays$RNA@meta.data <- cbind(seuratObj@assays$RNA@meta.data, to_append)
rm(to_append)

#combine cells into bulk samples
pseudo <- AggregateExpression(seuratObj, assays = "RNA", return.seurat = FALSE, group.by = group_by)
pseudo <- pseudo$RNA

#create a metadata dataframe.
col_data <- data.frame(bulk_group = colnames(pseudo)) %>%
  separate(col = bulk_group, into = group_by, sep = "_", remove = FALSE)

#PCA analysis####
#start with PCA to check how the pseudobulk samples cluster
pca <- apply(pseudo, 2, scale)
rownames(pca) <- rownames(pseudo)
pca <- t(pca)
pca <- prcomp(pca)

pca_data <- pca$x %>%
  as_tibble(rownames = "bulk_group") %>%
  left_join(col_data, by = "bulk_group") 

plot_list <- list()
for(group in group_by){
  plot_list[[group]] <- pca_data %>%
    ggplot(aes(x = PC1, y = PC2, color = .data[[group]], label = .data[[group]]))+
    geom_text()+
    ggtitle(group)
}

pca_plot <- ggalign::align_plots(!!!plot_list, nrow = 1)
pca_plot

#DE analysis####
##create DEseq object
dds <- DESeqDataSetFromMatrix(countData = as.data.frame(pseudo), colData = col_data, design = as.formula(paste0("~ ", group_by[1])))
#set the reference group
if(!is.na(reference_group)){
  dds[[group_by[1]]] <- relevel(dds[[group_by[1]]], ref = reference_group)
}
#run the analysis
dds <- DESeq(dds)

##Vulcano Plot####
plot_list <- list()
for(pair in resultsNames(dds)[2:length(resultsNames(dds))]){
  plot_list[[pair]] <- results(dds, name = pair) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    mutate(pvalue = padj) %>%
    mutate(color = pvalue < pvalue_threshold) %>%
    mutate(color = ifelse(color & log2FoldChange > log2FC_threshold, "blue", 
                          ifelse(color & log2FoldChange < -log2FC_threshold, "red", "#cccccc"))) %>%
    mutate(label = ifelse(color != "#cccccc", gene, NA)) %>%
    ggplot(aes(x = log2FoldChange, y = log10(pvalue), color = color, label = gene_names[label]))+
    geom_point()+
    geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), color = "#cccccc", linetype = "dashed")+
    geom_hline(yintercept = log10(pvalue_threshold), color = "#cccccc", linetype = "dashed")+
    geom_label_repel(size = 2)+
    guides(color = "none")+
    labs(title = pair)+
    scale_color_manual(values = c("red" = "red", "#cccccc" = "#cccccc", "blue" = "blue"))+
    scale_y_reverse()+
    theme_bw()+
    theme(panel.grid = element_blank())
}

vulcano_plot <- ggalign::align_plots(!!!plot_list, nrow = 1)
vulcano_plot

##Venn diagram####
#now lets check for overlapping DE genes 
#venn diagram
DE_list <- list()
for(pair in resultsNames(dds)[2:length(resultsNames(dds))]){
  res <- results(dds, name = pair) %>%
    as.data.frame() %>%
    filter(padj < pvalue_threshold) 
  
  DE_list$up[[pair]] <- res %>%
    filter(log2FoldChange > log2FC_threshold) %>%
    rownames()
  
  DE_list$down[[pair]] <- res %>%
    filter(log2FoldChange < -log2FC_threshold) %>%
    rownames()
}

venn_colors <- viridis::turbo(4)
up <- ggvenn::ggvenn(DE_list$up, fill_color = venn_colors)+ggtitle("Up-regulated")
down <- ggvenn::ggvenn(DE_list$down, fill_color = venn_colors)+ggtitle("Down-regulated")

venn_diagram <- ggalign::align_plots(up, down, nrow = 1)
venn_diagram