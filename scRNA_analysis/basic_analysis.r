#clean the environment####
rm(list = ls(all.names = TRUE))
gc()

#libraries####
library(Seurat)
library(tidyverse)
source("git_modules/scRNA_tools/scRNA_functions.R")
source("scRNA_analysis/env.R")

#set paths#### 
cellranger_outs_path <- "../inputs/scRNA/cellranger_outputs/sandfly_data/"
out_seurat_path <- "../inputs/scRNA/seurat_objects/sandfly_data/sandfly_processed.rds"

dirs <- list.dirs(cellranger_outs_path, recursive = FALSE)

#set parameters####
pca_dims <- 1:14

#read counts into seurat objects####
objects <- list()
for(dir in dirs){
  library <- basename(dir)
  counts <- Read10X(paste0(dir, "/filtered_feature_bc_matrix"), gene.column = 1)
  objects[[library]] <- CreateSeuratObject(counts, min.cells = 2, min.features = 200, project = library)
}

#remove outlier cells####
objects <- lapply(objects, function(x){
  x[["percent_mito"]] <- PercentageFeatureSet(x, pattern = "LdBPK-37")
  x <- subset(x, !percent_mito %in% boxplot.stats(percent_mito)$out) #removes outliers for percentage of mitochondrial RNA
  x <- subset(x, !nCount_RNA %in% boxplot.stats(nCount_RNA)$out) #removes outliers for percentage of mitochondrial RNA
  return(x)
})

#remove non-amastigotes from day 0####
#for that we map the day 0 library to the L.maj reference and keep only the cells flagged as amastigotes
#start by getting a L. donovani to L. major gene convertion table
Ld_to_Lmj <- xlsx::read.xlsx(Lmj_to_Ld, sheetIndex = 1) %>% 
  select(seurat_id, Lm) %>%
  deframe()

#map the object to the atlas
objects$`107706-001-001` <- meta_from_ref(
  query_obj = objects$`107706-001-001`,
  ref_obj = readRDS("../inputs/scRNA/seurat_objects/Lmj_Atlas.rds"),
  dims = pca_dims,
  metadata_to_extract = "stages.ident",
  prefix = "atlas.",
  id_conversion = Ld_to_Lmj,
  normalization_method = "SCT"
)

#remove non-amastigotes
objects$`107706-001-001` <- subset(objects$`107706-001-001`, atlas.stages.ident == "amastigotes")
objects$`107706-001-001`$atlas.stages.ident <- NULL
objects$`107706-001-001`$atlas.stages.ident.score <- NULL

#SCT transformation of each object####
objects <- lapply(objects, NormalizeData)
objects <- lapply(objects, ScaleData)

#merge the objects
seuratObj <- merge(objects[[1]], objects[(c(2:length(objects)))], add.cell.ids = names(objects), project = "sandflies")

#now run analysis on merged object
seuratObj <- JoinLayers(seuratObj)
seuratObj <- SCTransform(seuratObj, vars.to.regress = "percent_mito")
seuratObj <- RunPCA(seuratObj, verbose = FALSE)
seuratObj <- RunUMAP(seuratObj, dims = pca_dims)
seuratObj <- FindNeighbors(seuratObj, dims = pca_dims)
seuratObj <- FindClusters(seuratObj, resolution = 0.5)

#save the object to disk
saveRDS(seuratObj, file = out_seurat_path)
