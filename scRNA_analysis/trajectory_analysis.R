#trajectory inference script of the culture sample

##clean the environment####
rm(list = ls(envir=parent.frame()))
gc()

#libraries####
library(tidyverse)
library(cowplot)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(ggridges)
library(Seurat)
library(viridis)
library(slingshot)
library(tradeSeq)
library(RColorBrewer)
library(mclust)
library(hues)
library(parallel)
library(BiocParallel)
library(Cairo)
library(splines)
library(dtwclust)
library(clusterExperiment)
source("git_modules/scRNA_tools/scRNA_functions.R")

#set paths
seuratObj_path <- "../inputs/scRNA/seurat_objects/B12_E4_seurat_processed.rds"

#inputs####
seuratObj <- readRDS(seuratObj_path)
out_dir <- "../outputs/scRNA/trajectory_inference"

#set defaults (so these variables can be commented out below)
start_cluster <- NULL
end_cluster <- NULL
clusters_to_ignore <- NULL
seed_to_use <- NULL
lineages_to_test <- NULL

#set values
seed_to_use <- 3778473
traj_reduction <- "phate"
vis_reduction <- "phate"
dims_to_use <- 1:2
cluster_variable <- "atlas.stages.ident"
conditions_col <- "aneu_group"
start_cluster <- "procyclics"
#end_cluster <- "early_g1"
#end_cluster <- "g2_m"
#create a subdirectory where the trajectory results will be stored
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#to make results more reproducible, we set a seed and store it, so it can be reused later if needed.
if(is.numeric(seed_to_use)){
  seed_value <- seed_to_use
}else{
  warning("seting random seed...")
  seed_value <- sample(c(1:10000000), size = 1)
}
set.seed(seed_value)

#create GMM clusterings
if(cluster_variable == "GMM_clusters"){
  seuratObj@meta.data$GMM_clusters <- Mclust(Embeddings(seuratObj, reduction = traj_reduction, G = 1))$classification
}

#visualize the clustering
plot_number <- 1
plot <- DimPlot(seuratObj, reduction = traj_reduction, group.by = cluster_variable, pt.size = 0.25, label = TRUE, cols = seuratObj@misc$colors[[cluster_variable]])&
  ggtitle(paste("seed:", seed_value))
CairoPDF(paste0(out_dir, "/", plot_number,"_clusters_", toupper(traj_reduction), ".pdf"), width = 8, height = 6)
show(plot)
dev.off()
show(plot)
#delete any previous information on pseudotime from slingshot
seuratObj@meta.data[,grep("slingLineage_", colnames(seuratObj@meta.data))] <- NULL

#remove undesired clusters
cluster_data <- seuratObj@meta.data[,cluster_variable]
cluster_data[cluster_data %in% clusters_to_ignore] <- NA
##run slingshot on a dimentional reduction map.####
seuratObj@misc$slingshot <- slingshot(data = Embeddings(seuratObj, reduction = traj_reduction)[,dims_to_use], 
                                      clusterLabels = cluster_data,
                                      dist.method = "mnn",
                                      start.clus = start_cluster,
                                      end.clus = end_cluster, 
                                      shrink = 1,
                                      extend = "n",
                                      thresh = 0.01,
                                      stretch = 0.5,
                                      allow.breaks = FALSE)

#binding the Pseudotime values to the seurat metadata
slingData <- as.data.frame(slingPseudotime(seuratObj@misc$slingshot))
colnames(slingData) <- gsub("Lineage", "slingLineage_", colnames(slingData))
seuratObj@meta.data <- cbind(seuratObj@meta.data, slingData[rownames(seuratObj@meta.data),, drop = FALSE])
rm(slingData)

#plot all lineages found by slingshot
plot_list <- list()
for(lineage in names(slingLineages(seuratObj@misc$slingshot))){
  curvesData <- as.data.frame(slingCurves(seuratObj@misc$slingshot)[[lineage]]$s)
  colnames(curvesData) <- paste0("dim", c(1:ncol(curvesData))) 
  subtitle <- paste(seuratObj@misc$slingshot@metadata$lineages[[lineage]], collapse = " -> ")
  lineage <- gsub("Lineage", "slingLineage_", lineage)
  plot <- FeaturePlot(seuratObj, lineage,
                      reduction = traj_reduction)&
    geom_path(data = curvesData, aes(x = dim1, y = dim2), arrow = arrow(type = "closed"), linewidth = 1)&
    scale_color_gradientn(colors = rev(brewer.pal(11,'Spectral')[-6]), na.value = "#ededed")&
    ggtitle(lineage, subtitle)&
    theme(plot.title = element_text(hjust = 0), plot.subtitle = element_text(size = 6))
  plot_list[[lineage]] <- plot
}

#print lineage ordering
seuratObj@misc$slingshot@metadata$lineages

#combine the plots
plot_number <- plot_number + 1
plot <- plot_grid(plotlist = plot_list, nrow = 1)
CairoPDF(paste0(out_dir, "/", plot_number,"_lineages_", toupper(traj_reduction), ".pdf"), width = 2+(5*length(plot_list)), height = 5)
show(plot)
dev.off()
show(plot)


#set defaults
if(is.null(lineages_to_test)){
  lineages_to_test <- names(slingLineages(seuratObj@misc$slingshot))
  lineages_to_test <- as.integer(gsub("Lineage", "", lineages_to_test))
}

#plot just the selected lineages over the reduction
lineages <- paste0("Lineage", lineages_to_test)
curvesData <-slingCurves(seuratObj@misc$slingshot)[lineages] %>% 
  lapply(function(x)x$s) %>% 
  lapply(as.data.frame) %>%
  bind_rows(.id = "source")
  
plot <- DimPlot(seuratObj, reduction = traj_reduction, group.by = cluster_variable, cols = seuratObj@misc$colors[[cluster_variable]])
cols <- colnames(curvesData)[2:3]
for(lineage in unique(curvesData$source)){
  plot <- plot & geom_path(data = filter(curvesData, source == lineage), aes(x = .data[[cols[1]]], y = .data[[cols[2]]]), arrow = arrow(type = "closed"), linewidth = 1)
}

plot_number <- plot_number + 1
CairoPDF(paste0(out_dir, "/", plot_number,"_clusters_", toupper(traj_reduction), "_with_lineages.pdf"), width = 8, height = 6)
show(plot)
dev.off()
show(plot)
  
#check distribution of clusters along trajectories
to_plot <- seuratObj@meta.data %>%
  select(matches(c(conditions_col, cluster_variable)), contains("slingLineage")) %>%
  pivot_longer(cols = contains("slingLineage"), names_to = "lineage", values_to = "pseudotime") %>%
  mutate(lineage = gsub("slingLineage_", "Lineage ", lineage)) %>%
  filter(lineage %in% paste("Lineage", lineages_to_test))

#geom ridges remove NAs from the data before calculating densities, which artificially inflates distributions with lots of NAs.
#to fix it we set the NA values to double the maximum pseudotime in the data, and then we limit the plot scale to hide these values.
max_pseudotime <- max(to_plot$pseudotime, na.rm = TRUE)
plot <- to_plot %>%
  mutate(pseudotime = ifelse(is.na(pseudotime), max_pseudotime*2, pseudotime)) %>%
  ggplot(aes(x = pseudotime, 
             #y = fct_rev(.data[[cluster_variable]]), 
             #heigth = after_stat(density),
             fill = .data[[cluster_variable]]))+
  geom_density(alpha = 0.5, aes(y = after_stat(count)))+
  #geom_density_ridges(alpha = 0.5, panel_scaling = FALSE, scale = 5, from = 0, to = max_pseudotime, bandwidth = 0.3)+
  #scale_fill_manual(values = seuratObj@misc$colors[[cluster_variable]])+
  #labs(y = NULL, x = "Pseudotime", fill = "Annotation")+
  facet_grid(rows = vars(lineage), cols = vars(.data[[conditions_col]]))+
  #guides(fill = "none")+
  theme_ridges()+
  theme(panel.grid.major = element_line(color = "#f7f7f7"),
        text = element_text(size = 18))

plot_height <- length(unique(plot$data$lineage))*3+1
plot_widht <- 6
plot_number <- plot_number + 1
CairoPDF(paste0(out_dir, "/", plot_number,"_lineages_distributions.pdf"), width = plot_widht, height = plot_height)
show(plot)
dev.off()
show(plot)


#set parallelization backend
BPPARAM <- SnowParam(workers = 4, type = "SOCK")

#now evaluate k parameter
# evaluateK(counts = GetAssayData(seuratObj, layer = "counts", assay = "RNA"),
#           sds = seuratObj@misc$slingshot,
#           nGenes = 200, k = 3:12,
#           parallel = TRUE,
#           BPPARAM = BPPARAM)

## fit negative binomial GAMs to each gene####
seuratObj@misc$tradeSeq <- fitGAM(GetAssayData(seuratObj, layer = "counts", assay = "RNA"),
                               sds = seuratObj@misc$slingshot,
                               nknots = 6, 
                               parallel = TRUE,
                               BPPARAM = BPPARAM)

#save to disk
saveRDS(seuratObj, paste0(out_dir ,"/culture_sct_with_trajectory.rds"))

#create a named vector with gene names
gene_names <- setNames(seuratObj@assays$RNA@meta.data$cellranger_name, rownames(seuratObj@assays$RNA@meta.data))

#use the associationTest() function to find genes modulated along the pseudotime
association_test <- associationTest(seuratObj@misc$tradeSeq, global = FALSE, lineages = TRUE, l2fc = log2(2)) %>%
  as.data.frame()

# use the diffEndTest function to find DE genes at the end of the trajectories
diffEnd_test <- diffEndTest(seuratObj@misc$tradeSeq, global = FALSE, pairwise = TRUE, l2fc = log2(2)) %>%
  as.data.frame() 

#let's try the patternTest
pattern_test <- patternTest(seuratObj@misc$tradeSeq, global = FALSE, pairwise = TRUE, l2fc = log2(2)) %>%
  as.data.frame()
#for some reason pattern_test does not return the gene ids in the rownames. fix this here
rownames(pattern_test) <- rownames(seuratObj@misc$tradeSeq)

#also try startVsend
startVsEnd_test <- startVsEndTest(seuratObj@misc$tradeSeq, global = FALSE, lineages = TRUE, l2fc = log2(2)) %>%
  as.data.frame()

#colnames are different between different tradeseq test (really a consistency nightmare)
#fix this
colnames(startVsEnd_test) <- gsub("_lineage", "_", colnames(startVsEnd_test))
colnames(startVsEnd_test) <- gsub("FClineage", "FC_", colnames(startVsEnd_test))

#also use EarlyDEtest to identfy genes in the moment of differenciation 
earlyDE_test <- earlyDETest(seuratObj@misc$tradeSeq, global = FALSE, pairwise = TRUE, l2fc = log2(2), knots = c(4:5)) %>%
  as.data.frame()

#for some reason earlyDE_test also does not return the gene ids in the rownames. fix this here
rownames(earlyDE_test) <- rownames(seuratObj@misc$tradeSeq)

#create a list with the results
traj_results <- list(association_test = association_test, diffEnd_test = diffEnd_test, pattern_test = pattern_test, startVsEnd_test = startVsEnd_test, earlyDE_test = earlyDE_test)

#create vectors with the col names of the lineages
logFC_cols <- paste0("logFC_", lineages_to_test)
logFC_cols <- c(logFC_cols, paste0("logFC", apply(t(combn(lineages_to_test, 2)), 1, function(x)paste(x, collapse = "_"))))
pval_cols <- paste0("pvalue_", lineages_to_test)
pval_cols <- c(pval_cols, paste0("pvalue_", apply(t(combn(lineages_to_test, 2)), 1, function(x)paste(x, collapse = "vs"))))
wald_cols <- paste0("waldStat_", lineages_to_test)
wald_cols <- c(wald_cols, paste0("waldStat_", apply(t(combn(lineages_to_test, 2)), 1, function(x)paste(x, collapse = "vs"))))
df_cols <- paste0("df_", lineages_to_test)
df_cols <- c(df_cols, paste0("df_", apply(t(combn(lineages_to_test, 2)), 1, function(x)paste(x, collapse = "vs"))))
expected_cols <- c(wald_cols, df_cols, pval_cols, logFC_cols)
#adjust pvalues, select only the columns regarding the specified lineages, and set the order and filter parameters of each dataframe
traj_results <- lapply(traj_results, function(x){
  #remove unwanted columns
  x <- x[,colnames(x) %in% expected_cols]
  #adjust pvalues
  for(col in colnames(x)[colnames(x) %in% pval_cols]){
    x[,col] <- p.adjust(x[,col], method = "fdr")
  }
  #define the filter parameter (the minimum pvalue in the dataframe)
  x$filter_par <- apply(x[,colnames(x)[colnames(x) %in% pval_cols], drop = FALSE], 1, function(x){min(x, na.rm = TRUE)}) 
  
  #define the order parameter (either the wald or the log FC, whichever is present in the data)
  if(sum(grepl("logFC", colnames(x))) > 0){
    x$order_par <- apply(x[,colnames(x)[colnames(x) %in% logFC_cols], drop = FALSE], 1, function(x){max(abs(x), na.rm = TRUE)}) 
  }else{
    x$order_par <- apply(x[,colnames(x)[colnames(x) %in% wald_cols], drop = FALSE], 1, function(x){max(abs(x), na.rm = TRUE)}) 
  }
  
  #now we append the percentage of cells expressing those genes so we can filter low-expressed genes later
  x <- x %>%
    mutate(., mean_counts = rowMeans(counts(seuratObj@misc$tradeSeq)[rownames(.),])) %>%
    mutate(., percent_cells = rowMeans(counts(seuratObj@misc$tradeSeq)[rownames(.),] > 1))
  
  #now filter and arrange the datasets
  x <- x %>%
    filter(filter_par < 0.05) %>%
    filter(percent_cells > 0.01) %>%
    arrange(desc(order_par))
  return(x)
})

##plot heatmap with the genes####
#first get the pseudotime value for each cell in each lineage
pseudotime_data <- seuratObj@misc$slingshot@assays@data$pseudotime
#now for each gene set
for(i in c(1:length(traj_results))){
  result_name <- names(traj_results)[i]
  result_df <- traj_results[[i]]
  gene_set <- rownames(result_df)
  
  #get predicted smoother curves
  traj_data <- predictSmooth(seuratObj@misc$tradeSeq, gene = gene_set, nPoints = 100, tidy = TRUE) %>%
    #filter(lineage %in% lineages_to_test) %>%
    group_by(gene) %>%
    mutate(yhat = scale(yhat)) %>%
    ungroup()
  
  #cluster genes acording to temporal pattern expression
  gene_clusters <- traj_data %>%
    #filter(lineage == lineages_to_test[1]) %>%
    select(lineage, time, yhat, gene) %>%
    pivot_wider(names_from = c(lineage, time), values_from = yhat) %>%
    column_to_rownames("gene") %>%
    as.matrix() %>%
    proxy::dist(method = "dtw_basic") %>%
    hclust(method = "ward.D2")
  gene_groups <- cutree(gene_clusters, k = round(log10(length(gene_clusters$labels)))*2)
  
  #to order the genes we use the moment when they peak in expression along a smoother curve
  genes_order <- predictSmooth(seuratObj@misc$tradeSeq, gene = gene_set, nPoints = 10, tidy = TRUE) %>%
    group_by(gene) %>%
    filter(lineage == lineages_to_test[1]) %>%
    summarise(peak_expression = max(yhat, na.rm = T), peak_expression_time = time[which(yhat == peak_expression)]) %>%
    ungroup() %>%
    mutate(cluster = gene_groups[gene]) %>%
    arrange(peak_expression_time, cluster, peak_expression) %>%
    mutate(order = row_number()) %>%
    select(gene, order) %>%
    deframe()
  
  #append the gene group to the result data frame
  result_df <- result_df %>%
    rownames_to_column("gene_id") %>%
    mutate(gene_cluster = gene_groups[gene_id]) %>%
    mutate(gene_order = genes_order[gene_id])
  
  
  #export the dataframe as an excel
  colnames(result_df) <- gsub("_(\\d+)", "_lineage_\\1", colnames(result_df))
  
  result_df %>%
    select(-gene_order) %>%
    mutate(gene_name = seuratObj@assays$RNA@meta.data[gene_id,]$cellranger_name) %>%
    mutate(gene_id = gsub("-t42-1", "", gene_id)) %>%
    mutate(gene_id = gsub("-", "_", gene_id)) %>%
    mutate(TriTrypDB_link = paste0("https://tritrypdb.org/tritrypdb/app/record/gene/", gene_id)) %>%
    xlsx::write.xlsx(file = paste0(out_dir, "/", result_name, "_genes.xlsx"), row.names = FALSE)
  
  
  #create a color scale to represent the expression
  #set the range of values in the matrix. This is to ensure all heatmaps will have the same color scale
  color_breaks <- quantile(traj_data$yhat, c(0.01, 0.99))
  color_breaks <- seq(color_breaks[1], color_breaks[2], length = 11)
  color_vector <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  expression_colors <- colorRamp2(breaks = color_breaks, colors = color_vector)
  #create a color scale to represent the pseudotimes
  color_breaks <- sort(unique(traj_data$time))
  color_vector <- viridis::turbo(length(color_breaks))
  pseudotime_colors <- colorRamp2(breaks = color_breaks, colors = color_vector)
  
  #now for each lineage
  heatmap_list <- NULL
  for(lineage_to_plot in lineages_to_test){
    #first get the lineage data
    heat_data <- filter(traj_data, lineage == lineage_to_plot)
    #now we get the cells belonging to each pseudo time value, as well as their annotations
    cells_data <- pseudotime_data[,paste0("Lineage", lineage_to_plot)]
    cells_data <- data.frame(cell = names(cells_data), pseudotime = cells_data) %>%
      filter(!is.na(pseudotime)) %>%
      arrange(pseudotime) %>%
      mutate(morphotype = seuratObj@meta.data[cell,]$atlas.stages.ident)
    #now for each pseudotime interval in the smoothed data, we determine the average annotation of the cells in that interval
    pseudotime_intervals <- data.frame(pseudotime = unique(heat_data$time)) %>%
      mutate(next_time = lead(pseudotime)) 
    pseudotime_intervals$annotation <- NA
    for(j in c(1:(nrow(pseudotime_intervals)-1))){
     annotation <- cells_data %>%
       filter(pseudotime >= pseudotime_intervals$pseudotime[j]) %>%
       filter(pseudotime <= pseudotime_intervals$next_time[j]) %>%
       pull("morphotype") %>%
       as.character() %>%
       statistical_mode(na.rm = TRUE)
     pseudotime_intervals$annotation[j] <- annotation[1]
     rm(annotation)
    }
    pseudotime_intervals$annotation[nrow(pseudotime_intervals)] <- pseudotime_intervals$annotation[nrow(pseudotime_intervals)-1] 
    
    #now we convert the dataframe into a matrix for the heatmap 
    heat_data <- heat_data %>%
      select(gene, time, yhat) %>%
      pivot_wider(names_from = time, values_from = yhat) %>%
      column_to_rownames("gene") %>%
      as.matrix()
    
    #order the genes
    heat_data <- heat_data[names(genes_order),]
    
    #create the annotation of the pseudotime values with the corresponding morphotypes
    pseudotime_annotation <- as.numeric(colnames(heat_data))
    names(pseudotime_annotation) <- colnames(heat_data)
    morphotype_annotation <- pseudotime_intervals$annotation
    names(morphotype_annotation) <- pseudotime_intervals$pseudotime
    ann_col <- HeatmapAnnotation(pseudotime = pseudotime_annotation, morphotype = morphotype_annotation, col = list(morphotype = seuratObj@misc$colors$atlas.stages.ident, pseudotime = pseudotime_colors), show_annotation_name = lineage_to_plot == lineages_to_test[length(lineages_to_test)])
    heatmap <- Heatmap(heat_data, 
                       column_title = paste("Lineage", lineage_to_plot),
                       show_row_names = FALSE,
                       show_column_names = FALSE,
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       top_annotation = ann_col,
                       col = expression_colors,
                       use_raster = TRUE,
                       raster_quality = 10,
                       split = NULL,
                       heatmap_legend_param = list(title = "Scaled\nexpression"),
                       show_heatmap_legend = lineage_to_plot == lineages_to_test[1])
    heatmap_list <- heatmap_list + heatmap
  }
  #save the plot
  plot_number <- plot_number + 1
  plot_height <- 2+log2(nrow(heat_data))
  CairoPDF(paste0(out_dir, "/", plot_number,"_", result_name, "_trajectory_heatmap.pdf"),  width = 9, height = plot_height)
  show(heatmap_list)
  dev.off()
  
  #now plot the smoothers for each gene cluster
  #skip the smoother plot if there are too many genes to plot
  if(length(gene_set) > 600){
    next
  }
  

  for(cluster in unique(gene_groups)){
    #make the plot
    genes_to_plot <- names(gene_groups[which(gene_groups == cluster)])
    names(genes_to_plot) <- gene_names[genes_to_plot]
    genes_to_plot <- genes_to_plot[sort(names(genes_to_plot))]
    plot_width <- round(sqrt(length(genes_to_plot)))
    plot_height <- ceiling(length(genes_to_plot)/plot_width)
    plot_number <- plot_number + 1
    plot <- plot_smoothers(seuratObj@misc$tradeSeq, genes = genes_to_plot, lineages = lineages_to_test, gene_names = gene_names)&ggtitle(paste("Gene Cluster", cluster))
    CairoPDF(paste0(out_dir, "/", plot_number,"_", result_name, "_smoothers_GeneCluster_", cluster, ".pdf"), width = (plot_width*2)+4, height = (plot_height*1.5)+1)
    show(plot)
    dev.off()
  }
}

genes_to_plot <- rownames(seuratObj@assays$RNA@meta.data)[grep("meta domain|sherp|hydroph|surface antigen|gp63", tolower(seuratObj@assays$RNA@meta.data$cellranger_name))]
plot_smoothers(seuratObj@misc$tradeSeq, genes = sort(genes_to_plot), lineages = c(1,3), gene_names = gene_names, show_knots = TRUE)


