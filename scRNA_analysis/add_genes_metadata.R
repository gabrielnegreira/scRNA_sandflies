#this is a sub-script meant to be ran by another script.
#it will append relevant information to the metadata of the seurat object
source("scRNA_analysis/env.R") #paths to relevant information are set here

#get genome annotations
genome_annotation <- rtracklayer::readGFFAsGRanges(ref_annotation) %>%
  as.data.frame() %>%
  dplyr::rename(annotated_name = gene_name) %>%
  mutate(chromosome = as.numeric(factor(seqnames))) %>%
  mutate(seurat_id = gsub("_", "-", transcript_id))

#get TrypTag localizations
TrypTag_localizations <- read_delim(TrypTag_data) %>%
  column_to_rownames("Gene ID")

#get which genes are multicopy
multicopy_genes <- readxl::read_excel(ref_multicopy_genes) %>%
  mutate(seurat_id = gsub("_", "-", name))

#Tryptag has T.brucei genes. We need to replace them by L. donovani orthologs
TrypTag_localizations <- read_delim(Tb_to_Ld) %>%
  dplyr::rename(LdCL = source_id, Tb = `Input Ortholog(s)`) %>%
  select(LdCL, Tb) %>%
  separate_rows(Tb, sep = ",") %>% #create a 1-to-1 conversion table
  mutate(is_multicopy = LdCL %in% multicopy_genes$name) %>% #determine if a gene is multicopy in donovani
  mutate(LdCL = ifelse(is_multicopy, multicopy_genes$representative[match(LdCL, multicopy_genes$name)], LdCL)) %>% #replace multicopy genes by their representatives
  distinct() %>% #remove duplicates
  mutate(TrypTag_C_localization = TrypTag_localizations[Tb,]$`C localisation`) %>% #find the C localization
  mutate(TrypTag_N_localization = TrypTag_localizations[Tb,]$`N localisation`) %>% #find the N localization
  select(LdCL, TrypTag_C_localization, TrypTag_N_localization) %>%
  filter(!is.na(TrypTag_C_localization) | !is.na(TrypTag_N_localization)) %>% #remove entries with non-defined localization
  group_by(LdCL) %>%
  summarise(TrypTag_C_localization = paste(TrypTag_C_localization, collapse = ", "),
            TrypTag_N_localization = paste(TrypTag_N_localization, collapse = ", ")) %>%
  mutate(seurat_id = gsub("_", "-", LdCL))

#remove unwanted information
TrypTag_localizations <- TrypTag_localizations %>% select(-LdCL)

#combine the metadata into the seurat object
seuratObj[["RNA"]]@meta.data <- seuratObj[["RNA"]]@meta.data  %>%
  rownames_to_column("rowname") %>%
  left_join(genome_annotation, by = "seurat_id") %>%
  left_join(TrypTag_localizations, by = "seurat_id") %>%
  mutate(is_multicopy = seurat_id %in% multicopy_genes$seurat_id) %>%
  mutate(is_ribo = grepl("ribosom", tolower(annotated_name))) %>%
  mutate(is_mito = grepl("LdBPK-37", seurat_id)) %>%
  column_to_rownames("rowname")
