#control the environment for the scRNA analysis

#set paths (relative to repo root)####
##scripts####
cell_cycle_script <- "scRNA_analysis/score_cell_cycle.R"
flag_aneu_script <- "scRNA_analysis/flag_aneuploid_cells.R"


#reference genome annotations
ref_annotation <- "../inputs/reference_genome/LdCL_annotation_for_cellranger.gtf"
ref_multicopy_genes <- "../inputs/reference_genome/LdCL_multicopy_genes.xlsx"
ref_GO_annotations <- "../inputs/reference_genome/LdCL_GO_terms_updated.tsv"

#cross-species mapping
Tb_to_Ld <- "../inputs/external_data//Tb_to_LdCL_TriTrypDB.tsv"
Ltr_to_Ld <- "../inputs/external_data/LTR_to_LdCL.txt"

#orthologs for mapping to atlas
Lmj_to_Ld_sandflies <- "../inputs/reference_genome/Ld_sandfly_to_Lmj_one_to_one.tsv"
Lmj_to_Ld_culture <- "../inputs/reference_genome/Ld_culture_to_Lmj_one_to_one.tsv"


#data from other species
TrypTag_data <- "../inputs/external_data/TrypTag_localizations.tsv"
cell_cycle_tb_markers <- "../inputs/external_data/Tb_cell_cycle_genes.xlsx"
Lmj_ref_obj <- "../inputs/scRNA/seurat_objects/Lmj_Atlas.rds"
