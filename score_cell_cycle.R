source("custom_functions.R")
#source("2_set_common_variables.R")
##cell cycle scoring####
#link Tb to LdCl
expressed <- rowMeans(seuratObj)
expressed <- expressed[expressed > 0.1]
#get list of cell cycle genes from T.brucei with LdCL orthologs
tb_cell_cycle <- xlsx::read.xlsx("common_inputs/T_brucei_data/Tb_cell_cycle_genes.xlsx", sheetIndex = 1)
#get the list with the representatives of multicopy genes
multicopy_genes <- readxl::read_excel("common_inputs/LdCL_ref/LdCL_multicopy_genes.xlsx")

#create a list with genes for each cell cycle phase
cell_cycle_genes <- readr::read_delim("common_inputs/T_brucei_data/Tb_to_LdCL.tsv", delim = "\t") %>%
  dplyr::rename(LdCL = source_id, Tb = `Input Ortholog(s)`) %>%
  select(LdCL, Tb) %>%
  separate_rows(Tb, sep = ",") %>%
  filter(Tb %in% tb_cell_cycle$Tb_gene) %>%
  mutate(is_multicopy = LdCL %in% multicopy_genes$name) %>%
  mutate(LdCL = ifelse(is_multicopy, multicopy_genes$representative[match(LdCL, multicopy_genes$name)], LdCL)) %>% #replace multicopy genes by their representatives
  distinct() %>%
  mutate(cell_cycle_phase = tb_cell_cycle$cell_cycle_phase[match(Tb, tb_cell_cycle$Tb_gene)]) %>%
  mutate(LdCL = gsub("_", "-", LdCL)) %>%
  filter(LdCL %in% names(expressed)) %>%
  dplyr::select(cell_cycle_phase, LdCL) %>%
  split(.$cell_cycle_phase, .$LdCL) %>%
  sapply(function(x) x <- x$LdCL)

#assign cell cycle as in Brigs et al 2021
for(phase in names(cell_cycle_genes)){
  seuratObj <- MetaFeature(seuratObj, features = cell_cycle_genes[[phase]], meta.name = phase)
}
seuratObj@meta.data[,names(cell_cycle_genes)] <- seuratObj@meta.data[,names(cell_cycle_genes)]/colMeans(seuratObj@meta.data[,names(cell_cycle_genes)])
seuratObj@meta.data$non_cycling <- 1.5
cell_cycle_genes[["non_cycling"]] <- NA
seuratObj@meta.data$cell_cycle_phase <- names(cell_cycle_genes)[apply(seuratObj@meta.data[,names(cell_cycle_genes)], 1, which.max)]
#scCustomize::FeaturePlot_scCustom(seuratObj, features = names(cell_cycle_genes)[1:4], colors_use = umap_expr_colors)

seuratObj@meta.data <- seuratObj@meta.data %>%
  mutate(cell_cycle_phase = factor(cell_cycle_phase, levels = c("non_cycling", "early_g1", "late_g1", "s_phase", "g2_m")))

#create colors for cell cycle
seuratObj@misc$colors$cell_cycle_phase <- setNames(c("#bebebe", "#f8766d", "#7cae00", "#c77cff", "#01bfc4"),  c("non_cycling", "early_g1", "g2_m", "s_phase", "late_g1"))
