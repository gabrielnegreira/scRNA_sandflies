#this script does the general analysis, like dimensionality reduction, DE genes, etc.

#prepare the environment####
rm(list = ls(all.names = TRUE))

##load libraries####
library(glmGamPoi)
library(Seurat)
library(MAST)
library(tidyverse)
library(clustree)
library(Cairo)
library(ggalluvial)
library(ggrepel)
library(ggalign)
source("git_modules/scRNA_tools/scRNA_functions.R")
source("git_modules/GO_tools/GO_tools.R")

#get paths to external scripts and sources
source("scRNA_analysis/env.R")
Lmj_to_Ld <- Lmj_to_Ld
Lmj_ref_obj <- Lmj_ref_obj

#set paths
seuratObj_path <- "../inputs/scRNA/seurat_objects/seuratObj_all_timepoints_unprocessed.rds"
out_dir <- "../outputs/scRNA/scran_normalization_with_sct_no_integration"

##set parameters####
resolutions_to_try <- c(0:10)/10
chosen_resolution <- 0.3
pca_dims <- 1:20
norm_method <- "scran" #either scran or LogNormalize
use_sct <- TRUE #TRUE or FALSE
integrate_layers <- FALSE #TRUE or FALSE
categories_to_visualize <- c("pbm_day", "cell_cycle_phase", "seurat_clusters", "atlas.stages.ident")
qc_vars <- c("percent_mito", "nCount_RNA")
vars_to_regress <- c("percent_mito")
vars_to_regress <- NULL
n_variable_features <- 2000

#set outputs paths based on input object name####
base_name <- gsub("_unprocessed|.rds", "", basename(seuratObj_path))
out_seurat_name <- paste0(base_name, "_processed.rds")

##fix issues with memory limit####
library(future)
options(future.globals.maxSize = (32*(1024^3)))

##create output directory if needed####
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#actual analysis####
##read the seurat object####
seuratObj <- readRDS(file = seuratObj_path)

##remove outlier cells and lowly expressed genes####
####update cells metadata to include expression of mitochondrial and ribosomal genes
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

###split the object####
seuratObj[["RNA"]] <- split(seuratObj[["RNA"]], f = seuratObj$orig.ident)

###Normalize and scale genes####
if(norm_method == "scran"){
  seuratObj <- scran_normalize(seuratObj)
}else{
  seuratObj <- NormalizeData(seuratObj, normalization.method = "LogNormalize", scale.factor = 1e4)
}

if(use_sct){
  seuratObj <- SCTransform(seuratObj, vars.to.regress = vars_to_regress, variable.features.n = n_variable_features)
 }else{
  seuratObj <- FindVariableFeatures(seuratObj, nfeatures =  n_variable_features)
  seuratObj <- ScaleData(seuratObj) 
}

##run PCA####a
seuratObj <- RunPCA(seuratObj, verbose = FALSE)

### Integrate layers ####
if(integrate_layers){
  seuratObj <-  IntegrateLayers(
    object = seuratObj, 
    method = CCAIntegration, 
    normalization.method = c(SCT = "SCT", RNA = "LogNormalize")[DefaultAssay(seuratObj)])
  dr_to_use <- "integrated.dr"
}else{
  dr_to_use <- "pca"
}

##cluster cells####
set.seed(123)
seuratObj <- FindNeighbors(seuratObj, reduction = dr_to_use, dims = pca_dims, nn.method = "annoy") %>%
  FindClusters(resolution = resolutions_to_try)
seuratObj$seurat_clusters <- seuratObj[[paste0(DefaultAssay(seuratObj), "_snn_res.",chosen_resolution)]]
seuratObj <- SetIdent(seuratObj, value = "seurat_clusters")

##Run visualization DRs####
seuratObj <- RunUMAP(seuratObj, reduction = dr_to_use, dims = pca_dims)
#seuratObj <- RunPHATE(seuratObj)

#join layers####
seuratObj[["RNA"]] <- JoinLayers(seuratObj[["RNA"]])
seuratObj <- FindVariableFeatures(seuratObj, nfeatures =  n_variable_features) #has to be reran after joining

##define cell cycle####
source(cell_cycle_script)

##define morphotypes based on cell atlas####
if("atlas.stages.ident" %in% categories_to_visualize){
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
  seuratObj$atlas.stages.ident <- factor(seuratObj$atlas.stages.ident, 
                                         levels = c("procyclics", "procyclics-to-nectomonads", "nectomonads", "leptomonads", "haptomonads", "early-metacyclics", "late-metacyclics", "amastigotes"))
  seuratObj@misc[["colors"]]$atlas.stages.ident <- create_colors(seuratObj$atlas.stages.ident)
}

#Make plots####
##create color pallets for each category####
for(Category in categories_to_visualize){
  #get the category factor levels
  groups <- levels(as.factor(seuratObj@meta.data[[Category]]))
  
  #create a color scale if not present in misc
  if(
    is.null(seuratObj@misc[["colors"]][[Category]]) | 
    any(! names(seuratObj@misc[["colors"]][[Category]]) %in% groups)
  ){
    message(paste0("No colors found for category `", Category, "` in `misc` slot of the seurat object."))
    message("creating a new color scale for it using `create_colors`")
    seuratObj@misc[["colors"]][[Category]] <- create_colors(groups)
  }
}


##QC plots####
###qc variables####
qc_vars_plot <- seuratObj@meta.data %>%
  select(any_of(qc_vars), orig.ident) %>%
  pivot_longer(cols = -orig.ident, names_to = "feature") %>%
  group_by(orig.ident, feature) %>%
  mutate(value = value + rnorm(n(), 0, 1e-5)) %>% #this adds some noise to the values so violins look better
  ggplot(aes(x = orig.ident, y = value, group = orig.ident, fill = orig.ident))+
  geom_jitter(size = 0.03, show.legend = FALSE)+
  geom_violin(scale = "width", alpha = 0.5)+
  labs(x = NULL, y = NULL)+
  facet_wrap(vars(feature), scales = "free", nrow = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = paste0(out_dir, "/qc_vars_plot.pdf"),
  width = length(qc_vars) * 5 + 2,
  height = 6,
  scale = 1
)

###check if regressed vars still influence PCA####
always_plot <- c("nCount_RNA", "nFeature_RNA")
to_plot <- unique(c(always_plot, vars_to_regress))
for(regressed_var in to_plot){
  df_to_plot <- Embeddings(seuratObj, reduction = "pca")[,pca_dims] %>%
    as.data.frame() %>%
    rownames_to_column("cell_barcode") %>% 
    pivot_longer(cols = -cell_barcode, names_to = "component", values_to = "embedding") %>%
    left_join(
      seuratObj@meta.data %>%
        select(any_of(regressed_var)) %>%
        rownames_to_column("cell_barcode")
      , by = "cell_barcode"
    )%>%
    mutate(component = gsub("PC_", "", component)) %>%
    mutate(component = as.integer(component))
  
  #determine if the variable was regressed or not
  regressed_title <- ifelse(regressed_var %in% vars_to_regress, "(regressed)", "(not regressed)")
  
  plot <- ggplot(df_to_plot, aes(x = .data[[regressed_var]], y = embedding))+
    geom_point(size = 0.1)+
    geom_smooth(se = FALSE)+
    ggtitle(paste0("Effect of `", regressed_var, "` ", regressed_title," on PCA components"))+
    facet_wrap(vars(component))
  
  ggsave(filename = paste0(out_dir, "/", regressed_var, "VsPCA.pdf"), plot = plot, scale = 2)
}
 

### plot Variable Features ####
if(DefaultAssay(seuratObj) == "SCT"){
  variable_features_plot <- VariableFeaturePlot_SCT(seuratObj, label_info = "cellranger_name")
}else{
  variable_features_plot <- VariableFeaturePlot(seuratObj)
  variable_features_plot <- LabelPoints(plot = variable_features_plot, points =  head(VariableFeatures(seuratObj), 20), repel = TRUE)
}
ggsave(
  filename = paste0(out_dir, "/variable_features.pdf"),
  plot = variable_features_plot,
  scale = 2
)

### plot PCA Elbow####
elblow_plot <- ElbowPlot(seuratObj, ndims = max(pca_dims) + 20, reduction = "pca")+
  geom_vline(xintercept = max(pca_dims), linetype = "dashed", color = "darkred")+
  ggtitle(paste0("PCA Vectors from: ", "pca"),
          subtitle = paste0("assay: ", DefaultAssay(seuratObj)))

ggsave(
  filename = paste0(out_dir, "/elbow_plot.pdf"),
  plot = elblow_plot
)

### check cluster behaviour along different resolutions####
clustree_plot <- clustree(seuratObj)
clustree_plot <- clustree_plot+
  geom_hline(yintercept = length(resolutions_to_try) - which(resolutions_to_try == chosen_resolution), linetype = "dashed")+
  ggtitle(paste0("Chosen resolution: ", chosen_resolution, " (dashed line)"))

ggsave(
  filename = paste0(out_dir, "/clustree_plot.pdf"),
  plot = clustree_plot
)

##Main plots####
###Plot UMAPs grouped by categories
plot_list <- list()
for(Category in categories_to_visualize){
  #plot single UMAP grouped by category
  plot_list[[Category]] <- DimPlot(seuratObj, reduction = "umap", group.by = Category, cols = seuratObj@misc[["colors"]][[Category]])
}

ggsave(
  filename = paste0(out_dir, "/","umap_grouped_by_categories.pdf"),
  plot = ggalign::align_plots(!!!plot_list, ncol = 1),
  height = length(plot_list) * 4,
  width = 7,
  scale = 1.5
)

###Plot UMAPs splitted for each group category####
plot_list <- list()
for(Category in categories_to_visualize){
  #plot UMAP splitted by category, grouped by seurat clusters
  if(Category != "seurat_clusters"){
    plot <- DimPlot(seuratObj, reduction = "umap", group.by = "seurat_clusters", split.by = Category, cols = seuratObj@misc[["colors"]][["seurat_clusters"]])
    ggsave(
      filename = paste0(out_dir, "/UMAP_clusters_split_by_", Category, ".pdf"),
      plot = plot,
      width = length(unique(seuratObj@meta.data[[Category]])) * 4 + 1,
      height = 5,
      scale = 1.5
      )
  }
  
  #get the category factor levels
  groups <- levels(as.factor(seuratObj@meta.data[[Category]]))
  
  #now plot dimensionality reduction for each group
  plot_list[[Category]] <- list()
  for(Group in groups){
    
    #get cells to be highlighted
    cells <- seuratObj@meta.data %>%
      filter(.data[[Category]] == Group) %>%
      rownames()
   
    #make the plot for that group of cells
    plot <- DimPlot(seuratObj, reduction = "umap", cells.highlight = cells, pt.size = 0.05, sizes.highlight = 0.05)&
      scale_color_manual(values = c("lightgrey", as.character(seuratObj@misc[["colors"]][[Category]][[Group]])))&
      guides(color = "none")&
      ggtitle(Group)&
      theme_void()
    
    #save it to plot list
    plot_list[[Category]][[Group]] <- plot
  }
  
  #combine plots
  plot_list[[Category]] <- ggalign::align_plots(!!!plot_list[[Category]], nrow = 1) +
    ggalign::layout_title(title = Category)
  
  #save the plots 
  ggsave(filename = paste0(out_dir, "/", Category, "_umap.pdf"),
         plot = plot_list[[Category]],
         width = length(groups) * 5,
         height = 5,
         scale = 1,
         limitsize = FALSE)
}

###visualize the proportion of cells in each category####
for(Category in categories_to_visualize){
  plot <- seuratObj@meta.data %>%
    rename(ident = any_of(Category)) %>%
    select(orig.ident, ident)%>%
    group_by(orig.ident, ident)%>%
    summarise(ncells = n())%>%
    group_by(orig.ident)%>%
    mutate(percent = ncells/sum(ncells)) %>%
    mutate(percent_text = paste0(round(percent*100, 1), "%"))%>%
    ggplot(aes(x = orig.ident, alluvium = ident, stratum = ident, y = percent, fill = ident, label = percent_text))+
    geom_alluvium(alpha = 0.8, decreasing = FALSE)+geom_text_repel(stat = "stratum", decreasing = FALSE, direction = "y", nudge_x = 0.1)+
    #geom_col()+geom_text(position = position_stack(vjust = 0.5))+
    labs(fill = Category, x = "Library", y = "Percentage")+
    scale_fill_manual(values = seuratObj@misc[["colors"]][[Category]])+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(paste0(out_dir, "/", Category, "_proportions.pdf"), plot)
}

#plot flow of categories between each other
#plot category vs clusters
for(i in seq_along(categories_to_visualize)){
  for(j in seq_along(categories_to_visualize)){
    if(j >= i){
      next
    }
    cat1 <- categories_to_visualize[i]
    cat2 <- categories_to_visualize[j]
    
    plot <- seuratObj@meta.data %>%
      rename(ident1 = any_of(cat1), ident2 = any_of(cat2)) %>%
      select(ident1, ident2) %>%
      group_by(ident1, ident2) %>%
      summarise(n_cells = n()) %>%
      ggplot(aes(axis1 = ident1, axis2 = ident2, y = n_cells))+
      geom_alluvium(aes(fill = ident1))+
      geom_stratum(aes(fill = after_stat(stratum)))+
      geom_text(stat = "stratum", aes(label = after_stat(stratum)))+
      scale_fill_manual(values = c(seuratObj@misc[["colors"]][[cat1]], seuratObj@misc[["colors"]][[cat2]]))+
      scale_x_discrete(limits = c("seurat_clusters", "ident"),
                       labels = c(cat1, cat2),
                       expand = c(0.15, 0.05))+
      guides(fill = "none")+
      labs(y = "Number of cells")
    
    ggsave(paste0(out_dir,"/", cat1, "_Vs_", cat2, ".pdf"), plot)
  }
}

#check marker genes####
if(DefaultAssay(seuratObj) == "SCT"){
  seuratObj <- PrepSCTFindMarkers(seuratObj)
}
for(Category in categories_to_visualize){
  marker_genes <- FindAllMarkers(seuratObj, group.by = Category, logfc.threshold = log2(2), test.use = "MAST", only.pos = T, min.pct = 0.1)
  marker_genes <- marker_genes %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    as.data.frame() %>%
    mutate(transcript_id = gene) %>%
    dplyr::rename(gene_id = gene) %>%
    mutate(gene_id = gsub("LdCL-", "LdCL_", gene_id)) %>%
    mutate(gene_id = gsub("LdBPK-", "LdBPK_", gene_id)) %>%
    mutate(gene_id = gsub("-t42-1", "", gene_id)) %>%
    mutate(trypDB_link = ifelse(grepl("LdCL", gene_id), paste0("https://tritrypdb.org/tritrypdb/app/record/gene/", gene_id), NA)) %>%
    left_join(seuratObj[["RNA"]]@meta.data %>% select(-transcript_id) %>% rownames_to_column("transcript_id"), by = "transcript_id", keep =  FALSE)
  
  #save excel file with marker_genes
  xlsx::write.xlsx(marker_genes, paste0(out_dir, "/", Category, "_markers.xlsx"))
  
  genes_to_plot <- marker_genes %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 20) %>%
    pull("transcript_id") %>%
    unique()
  
  #plot heatmap####
  gene_names <- seuratObj@assays$RNA@meta.data %>%
    rownames_to_column("rowname") %>%
    select(rowname, cellranger_name) %>%
    deframe()
  
  plot <- DoHeatmap_plus(seuratObj, 
                         features = genes_to_plot,
                         feature_labels = gene_names,
                         group.by = c(Category, categories_to_visualize[categories_to_visualize != Category]),
                         group.colors = seuratObj@misc[["colors"]]
  )
  
  
  ggsave(plot = plot, filename = paste0(out_dir, "/", Category, "_markers_Heatmap.pdf"), height = 20, width = 30)
}

#save the seurat object to disk
saveRDS(seuratObj, file = paste0(dirname(seuratObj_path), "/", out_seurat_name))