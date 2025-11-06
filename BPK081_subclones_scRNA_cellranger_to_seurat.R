#this script is used to create the seurat Objects from the output files of cellranger-rna

#prepare the environment####
source("../bioinfo_tools/general_tools/clean_R_environment.R")

##load libraries####
library(Seurat)
library(rtracklayer)
library(tidyverse)
source("../bioinfo_tools/scRNA/scRNA_functions.R")

##set parameters####
project <- "BPK081_subclones"
cellranger_outputs <- "inputs/scRNA/cellranger_outputs/"
libraries <- c(`105118-001-001` = "B12", `105118-001-002` = "E4", `105118-001-003` = "G3", `105118-001-004` = "G6")
output_file <- "outputs/scRNA/BPK081_subclones_unprocessed.rds"

##fix issues with memory limit####
library(future)
options(future.globals.maxSize = (32*(1024^3)))

#code####
#create Seurat object####
##get the directory where data is
dirs <- list.dirs(cellranger_outputs, recursive = FALSE)
##create a list of Seurat objects
seuratObj_list <- list()
for(dir in dirs){
  library <- basename(dir)
  sample_name <- libraries[library]
  count_matrix <- Read10X(paste0(dir, "/filtered_feature_bc_matrix/"), gene.column = 1)
  seuratObj_list[[sample_name]] <- CreateSeuratObject(count_matrix, project = library)
}
##merge the list of seurat objects
seuratObj <- merge(seuratObj_list[[1]], seuratObj_list[2:4],  add.cell.ids = c(1:length(seuratObj_list)))
seuratObj <- JoinLayers(seuratObj)
rm(list = c("seuratObj_list", "count_matrix"))

##add the sample name to each cell
seuratObj@meta.data$subclone <- libraries[seuratObj@meta.data$orig.ident]

##create features metadata####
##get genome anotation
genome_annotation <- rtracklayer::readGFFAsGRanges("inputs/external_data/LdCL_ref/LdCL_annotation_for_cellranger.gtf") %>%
  as.data.frame() %>%
  dplyr::rename(annotated_name = gene_name) %>%
  mutate(chromosome = as.numeric(factor(seqnames)))

#get TrypTag localizations
TrypTag_localizations <- read_delim("inputs/external_data/T_brucei_data/TrypTag_localizations.tsv") %>%
column_to_rownames("Gene ID")
multicopy_genes <- readxl::read_excel("inputs/external_data/LdCL_ref/LdCL_multicopy_genes.xlsx")
TrypTag_localizations <- read_delim("inputs/external_data/T_brucei_data/Tb_to_LdCL.tsv") %>%
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
seuratObj[["RNA"]]@meta.data <- data.frame(cellranger_id = Read10X(paste0(cellranger_outputs,"/105118-001-001/filtered_feature_bc_matrix/"), gene.column = 1)@Dimnames[[1]],
                                           cellranger_name = Read10X(paste0(cellranger_outputs,"/105118-001-001/filtered_feature_bc_matrix/"), gene.column = 2)@Dimnames[[1]]) %>%
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

#tag aneuploid cells####
source("BPK081_subclones_scRNA_flag_aneuploid_cells.R")

#Map to Lmj atlas####
Ld_to_Lmj <- xlsx::read.xlsx("inputs/external_data/LdCL_ref/Lmj_to_LdCL_one_to_one.xlsx", sheetIndex = 1) %>%
  select(seurat_id, Lm) %>%
  deframe()
seuratObj <- meta_from_ref(query_obj = seuratObj, 
                           ref_obj = readRDS("inputs/external_data/reference_seurat_objects/Lmj_Atlas.rds"),
                           dims = 1:30,
                           subsample_ref = 5000,
                           metadata_to_extract = "stages.ident",
                           prefix = "atlas.",
                           normalization_method = "SCT",
                           ignore_values = list(stages.ident = "amastigotes"),
                           id_conversion = Ld_to_Lmj)

#create a color vector for the stages
seuratObj@misc$colors$atlas.stages.ident <- create_colors(seuratObj$atlas.stages.ident)

#store the seurat object to disk
saveRDS(seuratObj, file = output_file)
rm(seuratObj)

