#This script performs the basic seurat worflow and add morphotype anotations based on the L. major atlas.

#clean the environment####
rm(list = ls(all.names = TRUE))
gc()

#libraries####
library(Seurat)
library(tidyverse)
library(Matrix)
library(conflicted)
source("git_modules/scRNA_tools/scRNA_functions.R")
source("scRNA_analysis/env.R")
source("color_pallets.R")

#fix package conflicts
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rowSds", "MatrixGenerics")
#set paths#### 
cellranger_outs_path <- "../inputs/scRNA/cellranger_outputs/culture_data/"
out_seurat_path <- "../inputs/scRNA/seurat_objects/culture_data/culture_processed.rds"
out_plots_dir <- "../outputs/scRNA/culture_analysis/"

#get additional data####
Ld_to_Lmj <- read_tsv(Lmj_to_Ld_culture)

dirs <- list.dirs(cellranger_outs_path, recursive = FALSE)

#set parameters####
pca_dims <- 1:20
project_name <- "culture"

#read counts into seurat objects####
objects <- list()
for(dir in dirs){
  #get count matrix from cellranger outputs
  library <- basename(dir)
  counts <- Read10X(paste0(dir, "/filtered_feature_bc_matrix"), gene.column = 1)
  
  #create seurat object
  objects[[library]] <- CreateSeuratObject(counts, min.cells = 2, min.features = 200, project = library)
  
  #create metadata for features
  ##get the gene names set in cellranger (thus, from the genome annotation)
  cellranger_names <- Read10X(paste0(dir, "/filtered_feature_bc_matrix"), gene.column = 2) %>%
    rownames()
  
  #format it to include a list of ids compatible with seurat
  cellranger_names <- data.frame(cellranger_id = rownames(counts), cellranger_name = cellranger_names) %>%
    mutate(seurat_id = gsub("_", "-", cellranger_id)) %>%
    mutate(mean_counts = rowMeans(counts)[cellranger_id]) %>%
    mutate(sd_counts = MatrixGenerics::rowSds(counts)[cellranger_id])
  
  #create the metadata
  objects[[library]][["RNA"]]@meta.data <- data.frame(seurat_id = rownames(objects[[library]])) %>%
    left_join(cellranger_names, by = "seurat_id") %>%
    mutate(rowname = seurat_id) %>%
    column_to_rownames("rowname")
}

#remove outlier cells####
objects <- lapply(objects, function(x){
  x[["percent_mito"]] <- PercentageFeatureSet(x, pattern = "LdBPK-37")
  x <- subset(x, !percent_mito %in% boxplot.stats(percent_mito)$out) #removes outliers for percentage of mitochondrial RNA
  x <- subset(x, !nCount_RNA %in% boxplot.stats(nCount_RNA)$out) #removes outliers for percentage of mitochondrial RNA
  return(x)
})

#Normalize and scale each object####
#objects <- lapply(objects, NormalizeData)
#objects <- lapply(objects, ScaleData)

#merge the objects####
seuratObj <- objects[[1]]
seuratObj$project <- project_name

#now run analysis on merged object####
seuratObj <- SCTransform(seuratObj, vars.to.regress = "percent_mito", return.only.var.genes = FALSE)
seuratObj <- RunPCA(seuratObj, verbose = FALSE)
seuratObj <- RunUMAP(seuratObj, dims = pca_dims)
seuratObj <- FindNeighbors(seuratObj, dims = pca_dims)
seuratObj <- FindClusters(seuratObj, resolution = 0.6)

#assigned morphotypes to clusters/subclusters####
#start by getting a L. donovani to L. major gene convertion table
Ld_to_Lmj <- Ld_to_Lmj %>% 
  select(seurat_id, Lmj_ortholog) %>%
  distinct() %>%
  deframe()

##map it to reference
seuratObj <- meta_from_ref(
  query_obj = seuratObj,
  ref_obj = readRDS(Lmj_ref_obj),
  dims = pca_dims,
  metadata_to_extract = "stages.ident",
  prefix = "atlas.",
  id_conversion = Ld_to_Lmj
)

####get a table with the most frequent morphotype for each cluster
morpho_table <- seuratObj@meta.data %>%
  group_by(seurat_clusters, atlas.stages.ident) %>%
  summarise(n_cells = n()) %>%
  group_by(seurat_clusters) %>%
  mutate(n_cells_total = sum(n_cells)) %>%
  mutate(prop = n_cells/n_cells_total) %>%
  filter(prop == max(prop)) 

####since some clusters are mixture, we will find subclusters for them
#####first, define which clusters are mixed
mix_clusters <- morpho_table %>%
  filter(prop < 0.7) %>% #so any cluster where at least 30% of the cells were assigned a different morphotype is considered a mixed cluster
  pull("seurat_clusters") %>%
  as.character()

##### then for each mixed cluster, find subclusters
for(Cluster in mix_clusters){
  seuratObj <- FindSubCluster(seuratObj, cluster = Cluster, graph.name = paste0(DefaultAssay(seuratObj),"_snn"), subcluster.name = paste0("subcluster_of_", Cluster))
}

##### combine the results into a single column
cluster_cols <- grep("subcluster_of", colnames(seuratObj@meta.data), value = TRUE)
cluster_table <- seuratObj@meta.data[,cluster_cols, drop = FALSE]
seuratObj@meta.data$subcluster <- apply(cluster_table, 1, function(x){
  x <- as.character(x)
  counts <- table(x)
  
  if(length(counts) == 1){
    x <- x[1]
  }else{
    x <- names(counts)[counts == 1][1]
  }
  return(x)
})

####Now reassigned morphotypes
seuratObj@meta.data <- seuratObj@meta.data %>%
  rownames_to_column("rowname") %>%
  left_join(
    seuratObj@meta.data %>%
      group_by(seurat_clusters, subcluster, atlas.stages.ident) %>%
      summarise(n_cells = n()) %>%
      group_by(seurat_clusters) %>%
      #mutate(mother_cluster_morphotype = atlas.stages.ident[which(n_cells == max(n_cells))]) %>%
      group_by(seurat_clusters, subcluster) %>%
      filter(n_cells == max(n_cells)) %>%
      rename(assigned_morphotype = atlas.stages.ident) %>%
      group_by(subcluster) %>%
      filter(assigned_morphotype == assigned_morphotype[1]) %>% #in case of ties just chose one
      ungroup() %>%
      select(subcluster, assigned_morphotype) 
    , by = "subcluster"
  ) %>%
  column_to_rownames("rowname")

#set factor level for annotations####
seuratObj@meta.data$assigned_morphotype <- factor(seuratObj@meta.data$assigned_morphotype, levels = names(morphotype_colors))
seuratObj@meta.data$atlas.stages.ident <- factor(seuratObj@meta.data$atlas.stages.ident, levels = names(morphotype_colors))
seuratObj@meta.data$orig.ident <- factor(seuratObj@meta.data$orig.ident)

#define pbm_day####
seuratObj@meta.data <- seuratObj@meta.data %>%
  mutate(pbm_day = c("day00", "day02", "day04", "day06", "day10")[as.integer(orig.ident)])

#assign colors for clusters and morphotypes####
seuratObj@misc$colors$atlas.stages.ident <- morphotype_colors
seuratObj@misc$colors$assigned_morphotype <- morphotype_colors
seuratObj@misc$colors$pbm_day <- pbm_day_colors
seuratObj@misc$colors$seurat_clusters <- create_colors(unique(seuratObj@meta.data$seurat_clusters), palette = "all colors")
seuratObj@misc$colors$subcluster <- create_colors(unique(seuratObj$subcluster), palette = "fancy light")

##keep the colors for the unchanged clusters
unchanged_clusters <- which(names(seuratObj@misc$colors$subcluster) %in% names(seuratObj@misc$colors$seurat_clusters))
unchanged_clusters <- names(seuratObj@misc$colors$subcluster)[unchanged_clusters]
seuratObj@misc$colors$subcluster[unchanged_clusters] <- seuratObj@misc$colors$seurat_clusters[unchanged_clusters]

#append genes_metadata####

##first we need to fix feature metadata because it changes after JoinLayers
seuratObj[["RNA"]]@meta.data <- data.frame(seurat_id = rownames(seuratObj[["RNA"]])) %>%
  left_join(cellranger_names, by = "seurat_id") %>%
  mutate(mean_counts = rowMeans(seuratObj[["RNA"]]$counts)[seurat_id]) %>%
  mutate(sd_counts = rowSds(seuratObj[["RNA"]]$counts)[seurat_id]) %>%
  mutate(rowname = seurat_id) %>%
  column_to_rownames("rowname")

##now append extra information
source("scRNA_analysis/add_genes_metadata.R")

#visualize the annotations####
plot_list <- list()
for(Category in c("pbm_day", "atlas.stages.ident", "seurat_clusters", "assigned_morphotype")){
  plot_list[[Category]] <- DimPlot(seuratObj, group.by = Category, cols = seuratObj@misc$colors[[Category]])
}

final_plot <- ggalign::align_plots(!!!plot_list, ncol = 1)

dir.create(out_plots_dir, showWarnings = FALSE)
ggsave(paste0(out_plots_dir, "/UMAPs.pdf"), final_plot, height = 5 * length(plot_list), width = 8)

#save the object to disk####
saveRDS(seuratObj, file = out_seurat_path)
