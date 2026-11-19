#To be able to integrate our L. donovani dataset with the L. major atlas, I need a one-to-one conversion table between L. donovani genes and their L. major orthologs.
#This is challenging as some genes have no orthologs,while others have multiple. 
#This script exports a L.donovani to L.major conversion table using some heuristics to choose a single ortholog when multiple ones are present. 

#clean the environment####
rm(list = ls(all.names = TRUE))
gc()

#libraries#####
library(tidyverse)
library(xlsx)
library(Seurat)

#set path to store the conversion tables
out_dir <- "../inputs/reference_genome/"

#get LdCL orthologs (made with `make_ortholog_table.R`)
Ld_orthologs <- read_tsv("../inputs/reference_genome/TripDB_based_orthologs_for_LdCL.tsv")

#get the L. donovani count matrices
count_matrices <- list()
count_matrices[["Ld_sandfly"]] <- readRDS("../inputs/scRNA/seurat_objects/sandfly_data/sandfly_processed.rds")[["RNA"]]$counts
count_matrices[["Ld_culture"]] <- readRDS("../inputs/scRNA/seurat_objects/culture_data/culture_processed.rds")[["RNA"]]$counts

#get the L. major matrix
Lmj_matrix <- readRDS("../inputs/scRNA/seurat_objects/Lmj_Atlas.rds")[["RNA"]]@counts

#format the ortholog tables, one for each dataset
converstion_tables <- list()
for(dataset_name in grep("Ld", names(count_matrices), value = TRUE)){
  
  Ld_matrix <- count_matrices[[dataset_name]]
  
  converstion_table <- Ld_orthologs %>%
    filter(seurat_id %in% rownames(Ld_matrix)) %>% #removes L. donovani ids which are not present in the data
    mutate(Ld_mean_count = rowMeans(Ld_matrix)[seurat_id]) %>% #caculate mean counts of the Ld gene in the matrix (used for selection later)
    select(seurat_id, cellranger_id, Lmj_orthologs, Ld_mean_count) %>%
    separate_rows(Lmj_orthologs, sep = ",") %>%
    filter(Lmj_orthologs %in% rownames(Lmj_matrix)) %>% #removes any Lmj ortholog that by any change is not present in the atlas object
    mutate(Lmj_mean_count = rowMeans(Lmj_matrix)[Lmj_orthologs]) %>% #caculate mean counts of the Lmj ortholog in the reference matrix (used for selection later)
    group_by(cellranger_id) %>%
    mutate(n_Lmj_orthologs = n()) %>%
    mutate(to_keep = Lmj_mean_count == max(Lmj_mean_count, na.rm = TRUE)) %>% #if there are more than one ortholog, choose the one with the highers mean count in the atlas.
    filter(to_keep) %>%
    group_by(Lmj_orthologs) %>%
    mutate(n_Ld_orthologs = n()) %>%
    mutate(to_keep = Ld_mean_count == max(Ld_mean_count, na.rm = TRUE)) %>% #if multiple Ld genes map to the same Lmj gene, choose the Ld gene with highest mean count
    filter(to_keep) %>%
    group_by(cellranger_id) %>% 
    filter(Lmj_orthologs == Lmj_orthologs[1]) %>% #in the unlikely occasion of two homologs with the same average expression in both datasets is still present, just chose one
    group_by(Lmj_orthologs) %>% 
    filter(cellranger_id == cellranger_id[1]) %>% #in the unlikely occasion of two homologs with the same average expression in both datasets is still present, just chose one
    mutate(based_on_dataset = dataset_name) %>%
    rename(Lmj_ortholog = Lmj_orthologs) %>%
    select(seurat_id, cellranger_id, Lmj_ortholog, based_on_dataset)
  
  write_tsv(converstion_table, file = paste0(out_dir, "/", dataset_name,"to_Lmj_one_to_one.tsv"))
}