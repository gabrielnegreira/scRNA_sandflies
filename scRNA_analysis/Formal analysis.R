#this script does the standard analysis for the seurat object, including normalization, dimension reduction, clustering, etc

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

#solve package conflicts####
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("vars", "dplyr")
conflict_prefer("as.factor", "base")
conflict_prefer("align_plots", "ggalign")

#get paths to external scripts and sources
source("scRNA_analysis/env.R")
#variables set in `env.R`. Assignemnt here is just for code readability
Lmj_to_Ld <- Lmj_to_Ld 
Lmj_ref_obj <- Lmj_ref_obj

#get color pallets
source("color_pallets.R")
#just for readability:
morphotype_colors <- morphotype_colors
cell_cycle_colors <- cell_cycle_colors

#set paths
seuratObj_path <- "../inputs/scRNA/seurat_objects/seuratObj_all_timepoints_unprocessed.rds"
out_dir <- "../inputs/scRNA/seurat_objects/"

##set parameters####
qc_vars <- c("percent_mito", "nCount_RNA") #chose qc features to remove outlier cells
n_variable_features <- 850 #choose how many variable features to use (here, 10% of total genes)
pca_dims <- 1:13 #choose PCA dimensions to use for downstream analysis
cluster_resolution <- 0.3 #choose a resolution for cell clustering
subcluster_resolution <- 0.25 #choose a resolution for spiting clusters into subclusters

#set outputs paths based on input object name####
base_name <- gsub("_unprocessed|.rds", "", basename(seuratObj_path))
out_seurat_name <- paste0(base_name, "_processed.rds")

##fix issues with memory limit####
library(future)
options(future.globals.maxSize = (32*(1024^3))) #assumes 32 gb of ram available

##create output directory if needed####
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#actual analysis####
##read the seurat object####
seuratObj <- readRDS(file = seuratObj_path)

##remove outlier cells and lowly expressed genes####
###update cells metadata to include expression of mitochondrial and ribosomal genes
mito_features <- rownames(filter(seuratObj[["RNA"]]@meta.data, is_mitochondrial))
ribo_features <- rownames(filter(seuratObj[["RNA"]]@meta.data, is_ribosomal))

seuratObj <- PercentageFeatureSet(seuratObj, features = mito_features , col.name = "percent_mito")
seuratObj <- PercentageFeatureSet(seuratObj, features = ribo_features , col.name = "percent_ribo")

seuratObj@meta.data <- seuratObj@meta.data %>%
  rownames_to_column("rowname") %>%
  group_by(orig.ident) %>%
  mutate(
    across(
      all_of(qc_vars),
      ~ .x > quantile(.x, probs = 0.99, na.rm = TRUE),
      .names = "is_outlier_{.col}"
    )
  ) %>%
  ungroup() %>%
  mutate(to_remove = if_any(starts_with("is_outlier_"))) %>%
  column_to_rownames("rowname")

###remove outliers
seuratObj <- subset(seuratObj, !to_remove) 

##split the object per library####
seuratObj[["RNA"]] <- split(seuratObj[["RNA"]], f = seuratObj$orig.ident)

##Normalize and scale genes per library####
seuratObj <- scran_normalize(seuratObj)
seuratObj <- FindVariableFeatures(seuratObj, nfeatures =  n_variable_features)
seuratObj <- ScaleData(seuratObj) 

##run PCA####
seuratObj <- RunPCA(seuratObj, verbose = FALSE)

##cluster cells####
set.seed(123)
seuratObj <- FindNeighbors(seuratObj, reduction = "pca", dims = pca_dims, nn.method = "annoy")
seuratObj <- FindClusters(seuratObj, resolution = cluster_resolution)
seuratObj$seurat_clusters <- seuratObj[[paste0("RNA_snn_res.",cluster_resolution)]]
seuratObj <- SetIdent(seuratObj, value = "seurat_clusters")
seuratObj@misc$colors$seurat_clusters <- create_colors(unique(seuratObj$seurat_clusters))

##Run UMAP####
seuratObj <- RunUMAP(seuratObj, reduction = "pca", dims = pca_dims)

##join layers####
seuratObj[["RNA"]] <- JoinLayers(seuratObj[["RNA"]])
seuratObj <- FindVariableFeatures(seuratObj, nfeatures =  n_variable_features) #has to be reran after joining

##define cell cycle####
source(cell_cycle_script)
seuratObj@misc[["colors"]][["cell_cycle_phase"]] <- cell_cycle_colors

##map data to L major atlas to define morphotypes####
Ld_to_Lmj <- xlsx::read.xlsx(Lmj_to_Ld, sheetIndex = 1) %>% 
  select(seurat_id, Lm) %>%
  deframe()

seuratObj <- meta_from_ref(
  query_obj = seuratObj, 
  ref_obj = readRDS(Lmj_ref_obj),
  dims = pca_dims,
  #subsample_ref = 5000,
  metadata_to_extract = "stages.ident",
  prefix = "atlas.",
  normalization_method = ifelse(DefaultAssay(seuratObj) == "SCT", "SCT", "LogNormalize"),
  #vars_to_regress = "nCount_RNA",
  id_conversion = Ld_to_Lmj
)

#set factor levels
seuratObj$atlas.stages.ident <- factor(seuratObj$atlas.stages.ident, 
                                         levels = names(morphotype_colors))
#set color scale for it
seuratObj@misc$colors$atlas.stages.ident <- morphotype_colors

###assign morphotypes for each cluster####
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
  seuratObj <- FindSubCluster(seuratObj, cluster = Cluster, graph.name = "RNA_snn", resolution = subcluster_resolution, subcluster.name = paste0("subcluster_of_", Cluster))
}

##### combine the results into a single column
cluster_cols <- grep("subcluster_of", colnames(seuratObj@meta.data), value = TRUE)
seuratObj@meta.data$subcluster <- apply(seuratObj@meta.data[,cluster_cols, drop = FALSE], 1, function(x){
  x <- as.character(x)
  counts <- table(x)
  
  if(length(counts) == 1){
    x <- x[1]
  }else{
    x <- names(counts)[counts == 1][1]
  }
  return(x)
})

#####create color pallete for additional subclusters
seuratObj@misc$colors$subcluster <- create_colors(unique(seuratObj$subcluster), palette = "fancy light")
######keep the colors for the unchanged clusters
unchanged_clusters <- which(names(seuratObj@misc$colors$subcluster) %in% names(seuratObj@misc$colors$seurat_clusters))
unchanged_clusters <- names(seuratObj@misc$colors$subcluster)[unchanged_clusters]
seuratObj@misc$colors$subcluster[unchanged_clusters] <- seuratObj@misc$colors$seurat_clusters[unchanged_clusters]

####Now reassigned morphotypes
seuratObj@meta.data <- seuratObj@meta.data %>%
  rownames_to_column("rowname") %>%
  left_join(
    seuratObj@meta.data %>%
      group_by(subcluster, atlas.stages.ident) %>%
      summarise(n_cells = n()) %>%
      group_by(subcluster) %>%
      filter(n_cells == max(n_cells)) %>%
      rename(assigned_morphotype_final = atlas.stages.ident) %>%
      group_by(subcluster) %>%
      mutate(n = n()) %>% 
      filter(n == 1 | n > 1 & assigned_morphotype_final != "leptomonads") %>% #in case of a tie will favor other forms than leptomonads
      select(subcluster, assigned_morphotype_final) 
    , by = "subcluster"
  ) %>%
  column_to_rownames("rowname")

####set colors for new assigned morphotypes
seuratObj@misc$colors$assigned_morphotype_final <- seuratObj@misc$colors$atlas.stages.ident

#save the seurat object to disk
saveRDS(seuratObj, file = paste0(out_dir, "/", out_seurat_name))