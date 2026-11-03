#this script is used to create the seurat Objects from the output files of cellranger-rna
# it will also append several gene metadata based on external sources such as TrypTag,
#clean the environment####
rm(list = ls(all.names = TRUE))
gc()

##load libraries####
library(Seurat)
library(rtracklayer)
library(tidyverse)
source("git_modules/scRNA_tools/scRNA_functions.R")

#get required paths
source("scRNA_analysis/env.R")

##set parameters####
libraries <- c(
  `107706-001-001` = "day00", 
  `107706-001-002` = "day02", 
  `107706-001-003` = "day04", 
  `107706-001-004` = "day06",
  `107706-001-005` = "day10"
  )

## set paths####
cellranger_out_dir <- "../inputs/scRNA/cellranger_outputs/"
seurat_out_dir <- "../inputs/scRNA/seurat_objects/"

##fix issues with memory limit####
library(future)
options(future.globals.maxSize = (32*(1024^3)))

#code####
#create Seurat object####
##get the directory where data is
dirs <- list.dirs(cellranger_out_dir, recursive = FALSE)
#create a list of Seurat objects
seuratObj_list <- list()
for(dir in dirs){
  library <- basename(dir)
  sample_name <- libraries[library]
  count_matrix <- Read10X(paste0(dir, "/filtered_feature_bc_matrix/"), gene.column = 1)
  seuratObj_list[[sample_name]] <- CreateSeuratObject(count_matrix, project = library)
}
##merge the list of seurat objects
seuratObj <- merge(seuratObj_list[[1]], seuratObj_list[2:5],  add.cell.ids = c(1:length(seuratObj_list)))
seuratObj <- JoinLayers(seuratObj)
rm(list = c("seuratObj_list", "count_matrix"))

##add the sample name to each cell
seuratObj@meta.data$pbm_day <- libraries[seuratObj@meta.data$orig.ident]

##create features metadata####
##get genome anotation
genome_annotation <- rtracklayer::readGFFAsGRanges(ref_annotation) %>%
  as.data.frame() %>%
  dplyr::rename(annotated_name = gene_name) %>%
  mutate(chromosome = as.numeric(factor(seqnames)))

#get TrypTag localizations
TrypTag_localizations <- read_delim(TrypTag_data) %>%
column_to_rownames("Gene ID")
multicopy_genes <- readxl::read_excel(ref_multicopy_genes)
TrypTag_localizations <- read_delim(Tb_to_Ld) %>%
  dplyr::rename(LdCL = source_id, Tb = `Input Ortholog(s)`) %>%
  select(LdCL, Tb) %>%
  separate_rows(Tb, sep = ",") %>%
  mutate(is_multicopy = LdCL %in% multicopy_genes$name) %>%
  mutate(LdCL = ifelse(is_multicopy, multicopy_genes$representative[match(LdCL, multicopy_genes$name)], LdCL)) %>% #replace multicopy genes by their representatives
  distinct() %>%
  mutate(LdCL = gsub("_", "-", LdCL)) %>%
  mutate(TrypTag_C_localization = TrypTag_localizations[Tb,]$`C localisation`) %>%
  mutate(TrypTag_N_localization = TrypTag_localizations[Tb,]$`N localisation`) %>%
  select(LdCL, TrypTag_C_localization, TrypTag_N_localization) %>%
  filter(!is.na(TrypTag_C_localization) | !is.na(TrypTag_N_localization)) %>%
  group_by(LdCL) %>%
  summarise(TrypTag_C_localization = paste(TrypTag_C_localization, collapse = ", "),
            TrypTag_N_localization = paste(TrypTag_N_localization, collapse = ", ")) %>%
  column_to_rownames("LdCL")
  
#combine the metadata into the seurat object
seuratObj[["RNA"]]@meta.data <- data.frame(cellranger_id = Read10X(paste0(dirs[[1]], "/filtered_feature_bc_matrix/"), gene.column = 1)@Dimnames[[1]],
                                           cellranger_name = Read10X(paste0(dirs[[1]], "/filtered_feature_bc_matrix/"), gene.column = 2)@Dimnames[[1]]) %>%
  cbind(., genome_annotation[match(.$cellranger_id, genome_annotation$transcript_id),]) %>%
  mutate(seurat_id = gsub("_", "-", cellranger_id)) %>%
  mutate(is_mitochondrial = grepl("ldbpk_37", tolower(cellranger_id))) %>%
  mutate(is_ribosomal = grepl("ribosom", tolower(cellranger_name))) %>%
  mutate(is_multicopy = cellranger_id %in% multicopy_genes$name) %>%
  mutate(mean_count = rowMeans(seuratObj[["RNA"]]$counts, na.rm = TRUE)[seurat_id]) %>%
  mutate(sd_count = apply(seuratObj[["RNA"]]$counts, 1, sd, na.rm = TRUE)[seurat_id]) %>%
  mutate(ncells_express = rowSums(seuratObj[["RNA"]]$counts > 0, na.rm = TRUE)[seurat_id]) %>%
  mutate(prop_cells_express = rowMeans(seuratObj[["RNA"]]$counts > 0, na.rm = TRUE)[seurat_id]) %>%
  mutate(TrypTag_C_localization = TrypTag_localizations[seurat_id,]$TrypTag_C_localization) %>%
  mutate(TrypTag_N_localization = TrypTag_localizations[seurat_id,]$TrypTag_N_localization) %>%
  rownames_to_column("remove") %>%
  column_to_rownames("seurat_id") %>%
  select(-remove)

rm(list = c("dirs", "dir", "library", "libraries", "sample_name", "genome_annotation"))

#Map to Lmj atlas####
Ld_to_Lmj <- xlsx::read.xlsx(Lmj_to_Ld, sheetIndex = 1) %>%
  select(seurat_id, Lm) %>%
  deframe()
seuratObj <- meta_from_ref(query_obj = seuratObj, 
                           ref_obj = readRDS(Lmj_ref_obj),
                           dims = 1:30,
                           subsample_ref = 5000,
                           metadata_to_extract = "stages.ident",
                           prefix = "atlas.",
                           normalization_method = "SCT",
                           #ignore_values = list(stages.ident = "amastigotes"),
                           id_conversion = Ld_to_Lmj)

#create a color vector for the stages
seuratObj@misc$colors$atlas.stages.ident <- create_colors(seuratObj$atlas.stages.ident)

#write the object to disk
saveRDS(seuratObj, "../inputs/scRNA/seurat_objects/seuratObj_unprocessed.rds")