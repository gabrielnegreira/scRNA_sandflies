#This script exports a dataframe with the homologs defined for LdCL transcripts
#clean the environment
rm(list = ls(all.names = TRUE))
gc()

#libraries#####
library(tidyverse)
library(xlsx)
library(Seurat)

#get homolog tables exported from trypDB
ortho_tables <- list()
ortho_tables[["Lmj"]] <- read_tsv("../inputs/external_data//LmjF_to_LdCL_TriTrypDB.tsv")
ortho_tables[["Tb"]] <-  read_tsv("../inputs/external_data//Tb_to_LdCL_TriTrypDB.tsv") 
ortho_tables[["Ltr"]] <- read_tsv("../inputs/external_data//LTR_to_LdCL_TriTrypDB.tsv")

#format the orthology tables
for(species in names(ortho_tables)){
  ortho_col <- paste0(species, "_orthologs")
  
  ortho_tables[[species]] <- ortho_tables[[species]] %>%
    rename(cellranger_id = source_id, !!ortho_col := `Input Ortholog(s)`) %>%
    select(cellranger_id, all_of(ortho_col))
}
  
#add sherp orthologs
ortho_tables[["Lmj"]] <- rbind(ortho_tables[["Lmj"]], c("LdCL_230019430-t42_1", "LmjF.23.1050"))

#get the gene ids of the reference genome used in cellranger
Ld_ids <- Read10X("../inputs/scRNA/cellranger_outputs/sandfly_data/107706-001-001/filtered_feature_bc_matrix/", 1) %>%
  rownames() %>%
  data.frame(cellranger_id = .) %>%
  mutate(seurat_id = gsub("_", "-", cellranger_id)) #seurat converts underscores (_) to hifens (-). Fix it here.

for(species in names(ortho_tables)){
  ortho_table <- ortho_tables[[species]]
  Ld_ids <- Ld_ids %>%
    left_join(ortho_table, by = "cellranger_id")
  rm(ortho_table)
}

#export it 
write_tsv(Ld_ids, "../inputs/reference_genome/TripDB_based_orthologs_for_LdCL.tsv")