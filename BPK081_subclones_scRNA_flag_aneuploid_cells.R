#this script will tag cells as aneuploid or not based on differences bwtween the average expression of their chromosomes.
#It was designed to work with seurat objects as inputs. 
#In case prior knowledge on which chromosomes should be aneuploid in each sample is available, this can be provided as the aneu_chromo list. 
#in the feature metadata in the seurat object there should be a column named 'chromosome' indicating in which chromosome each gene is located, 
#as well as a 'start' and 'end' column indicating where each gene is located in the chromosome and a 'strand' column indicating in which stran (+ or -) each gene is located.

#libraries####
library(tidyverse)
library(ggplot2)
library(ggridges)
library(Seurat)
library(viridis)
library(Matrix)
library(factoextra)

#inputs####
#setwd("/Users/gnegreira/Dropbox/work/ITM/PhD/Scripts/Mine/R/scRNA data analysis/")
#seuratObj <- readRDS("outputs/ran_with_LdCL/BPK081_subclones/BPK081_subclones_seuratObj.rds")

#parameters####
cap_type <- "3" #is single-cell RNA sequencing done at 3' or 5'?
subsample = FALSE #If true, will subsample 10% of the cells to perform the analysis. Usefull to speed up the script while developing it. Not working currently
aneu_chromo <- list("B12" = c(1, 5, 10, 13, 20, 23),
                    "E4" = c(13, 20, 23, 26),
                    "G3" = c(20, 23, 26, 29),
                    "G6" = c(1, 9, 13, 16, 20, 22, 23, 24, 26, 35))
#aneu_chromo <- NULL
slot_to_use <- "counts"
group_cells_by <- "subclone"
ignore_variable_genes <- FALSE #will ignore variable genes if TRUE. Variable genes are defined by the FindVariableFeatures() function of Seurat 
chromosomes_to_ignore <- c(37) #a vector containing chromosomes to ignore. 
make_plots <- TRUE # if TRUE will plot several graphs along the analysis. This part will be improved later.
tag_method <- "aneu_ratio" #can be aneu_ratio or "pca". Only aneu_ratio is working for now. 

#custom_functions####
#this function will calculate the gini index of a numeric vector. I tried to use it to distinguish aneuploid cells but was not very helpful.
gini <- function (x, weights = rep(1, length = length(x)), na.rm = FALSE){
  x <- unlist(x)
  
  if (na.rm) {
    k <- is.na(x)
    if (any(k)) 
      x <- x[!k]
  }
  
  ox <- order(x)
  x <- x[ox]
  weights <- weights[ox]/sum(weights)
  p <- cumsum(weights)
  nu <- cumsum(weights * x)
  n <- length(nu)
  nu <- nu/nu[n]
  sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}

#code####
#get the metadata for the genes. 
genes_meta <- as.data.frame(seuratObj@assays$RNA@meta.data)
if(cap_type == "5"){
  warning("assuming 5' capture was used...")
  genes_meta <- genes_meta %>%
    mutate(cap_position = ifelse(strand == "+", start, end))
}else{
  warning("assuming 3' capture was used...")
  genes_meta <- genes_meta %>%
    mutate(cap_position = ifelse(strand == "+", end, start))
}

#get the metadata for the cells
cells_meta <- as.data.frame(seuratObj@meta.data)

#get the expression matrix. Should be a sparse matrix. 
expr_table <- seuratObj[["RNA"]][slot_to_use]

#removing NA chromosomes
genes_meta <- filter(genes_meta, !is.na(chromosome))
expr_table <- expr_table[rownames(genes_meta),]

#removing chromosomes to be ignored
if(length(chromosomes_to_ignore > 0)){
  genes_meta <- filter(genes_meta, !chromosome %in% chromosomes_to_ignore)
  expr_table <- expr_table[rownames(genes_meta),]
}

#removing variable genes
if(ignore_variable_genes){
  print("variable genes will be excluded from this analysis")
  genes_meta <- genes_meta[which(!rownames(genes_meta) %in% VariableFeatures(seuratObj)),]
  expr_table <- expr_table[rownames(genes_meta),]
}

#calculating average expression of each chromosome in each cell####
#get a list with all chromosomes
chromo_list <- genes_meta %>%
  filter(!is.na(chromosome)) %>%
  filter(!chromosome %in% chromosomes_to_ignore) %>%
  select(chromosome)%>%
  unique() %>%
  unlist() %>%
  as.character()

#create a new matrix where the average expression values will be stored
avg_chromo_matrix <- matrix(ncol = ncol(expr_table), nrow = length(chromo_list))
rownames(avg_chromo_matrix) <- chromo_list
colnames(avg_chromo_matrix) <- colnames(expr_table)

#for each chromosome
for(i in c(1:length(chromo_list))){
  print(i)
  #select the values of the genes of those chromosome
  values <- expr_table[as.character(rownames(filter(genes_meta, chromosome == chromo_list[i]))),]
  #remove outlier genes
  values[which(values > quantile(values, 0.999))] <- NA
  # values <- apply(values, 2, function(x){
  #   x[which(x > quantile(x, 0.99))] <- NA
  #   return(x)
  # })
  #calculate the means of values and append to the new matrix.
  avg_chromo_matrix[i,] <- colMeans(values, na.rm = T)
}


#normalize and scale the avg_chromo_matrix
avg_chromo_matrix <- t(t(avg_chromo_matrix)/colMeans(avg_chromo_matrix))
avg_chromo_matrix <- avg_chromo_matrix *2

#if a chromosome has no gene being expressed in a cell, this will result in a 0/0 operation, leading to an NA. Thus, we replace here NAs by 0.
avg_chromo_matrix[is.na(avg_chromo_matrix)] <- 0

#calculate the chromosomal variance and gini of each cell####
#transform the avg matrix into a long dataframe for downstream analysis
cells_aneu_meta <- avg_chromo_matrix %>%
  as.data.frame()%>%
  rownames_to_column("chromosome")%>%
  pivot_longer(cols = -"chromosome", names_to = "cell", values_to = "mean_expr") %>%
  group_by(cell) %>%
  summarize(karyo_gini = gini(mean_expr, na.rm = T),
            karyo_sd = sd(mean_expr, na.rm = T)) %>%
  column_to_rownames("cell")

#remove any homonimous column
cells_meta <- cells_meta[,which(!colnames(cells_meta) %in% colnames(cells_aneu_meta))]

#append the aneu_meta to the cells metadata
cells_meta <- cbind(cells_meta, cells_aneu_meta[match(rownames(cells_meta), rownames(cells_aneu_meta)),])

#calculate the ratio between polysomic and disomic chromosomes
#create a new column in the metadata where the aneu_ratio will be stored.
cells_meta$aneu_ratio <- NA
#get the sample names
samples <- unique(cells_meta[,group_cells_by])
#get the chromosome names
chromosomes <- rownames(avg_chromo_matrix)

plots <- list()
#for each sample
for(i in c(1:length(samples))){
  #define the sample
  sample_to_select <- samples[i]
  #define the cells from that sample
  cells_to_select <- cells_meta %>%
    filter(!!sym(group_cells_by) == sample_to_select)%>%
    rownames()
  #define the chromosomes which are expected to be aneuploid in that sample
  if(is.list(aneu_chromo)){
    #in case an aneu_chromo list is provided, use it to calculate the ratio between the average expression of disomic and polysomic chromosomes
    polysomic_chromosomes <- chromosomes[which(chromosomes %in% aneu_chromo[[sample_to_select]])]
    disomic_chromosomes <- chromosomes[which(!chromosomes %in% aneu_chromo[[sample_to_select]])]
  }else{
    #otherwise try to figure out which chromosomes are aneuploid in each sample
    warning("no aneu_chromo list provided. Will try to differentiate polysomic from disomic chromosomes using PCA and kmeans...")
    #get the values of the cells for that sample
    sample_matrix <- avg_chromo_matrix[, cells_to_select]
    #compress the data in a PCA
    pca <- prcomp(sample_matrix)
    #use kmeans to separate the chromosomes into two groups: still experimental. not sure how it will work with homogeneous data.
    chromo_clusters <- kmeans(scale(pca$x[,1:2]), centers = 2, nstart = 50)$cluster
    #create a dataframe with the summary for all chromosomes
    chromo_info <- data.frame(chromo = names(chromo_clusters), 
                                 cluster = chromo_clusters, 
                                 mean_expr =  rowMeans(sample_matrix, na.rm = TRUE))
    #now we need to know which cluster represents the disomic and which represent the polysomic ones.
    chromo_summary <- chromo_info %>%
      group_by(cluster) %>%
      summarise(mean = mean(mean_expr, na.rm = T)) %>%
      arrange(desc(mean)) %>%
      mutate(group = ifelse(mean == max(mean, na.rm = T), "polysomic", "disomic"))
    #the chromosome group with the highest mean expression will be assigned as 'polysomic'. The others will be disomic.
    chromo_info$group <- chromo_summary$group[match(chromo_info$cluster, chromo_summary$cluster)]
    polysomic_chromosomes <- chromo_info$chromo[which(chromo_info$group == "polysomic")]
    disomic_chromosomes <- chromo_info$chromo[which(chromo_info$group == "disomic")]
    #if polysomic_chromosomes has just one chromosome, it the colMeans function will not work. This is a workaround:
    if(length(polysomic_chromosomes) < 2){
      warning("found only one polysomic chromosome. Analysis is likely unreliable.")
      polysomic_chromosomes <- rep(polysomic_chromosomes, 2)
    }
    #same for disomic ones
    if(length(disomic_chromosomes) < 2){
      warning("found only one disomic chromosome. Analysis is likely unreliable.")
      disomic_chromosomes <- rep(disomic_chromosomes, 2)
    }
    if(make_plots){
      plots[[sample_to_select]] <- pca$x %>%
        as.data.frame()%>%
        rownames_to_column("chromo")%>%
        mutate(cluster = chromo_clusters[chromo]) %>%
        ggplot(aes(x = PC1, y = PC2, label = chromo, fill = factor(cluster)))+
        geom_label()+
        labs(title = paste("PCA and k-means clustering of chromosomes in sample:", sample_to_select),
             fill = "Cluster")
    }
  }
  
  #calculate the average expression of the polysomic chromosomes
  print(paste("using chromosomes", paste(polysomic_chromosomes, collapse = ", "), "to detect aneuploid cells in sample", sample_to_select))
  values_polysomic <- colMeans(avg_chromo_matrix[polysomic_chromosomes, cells_to_select])
  #calculate the average expression of the disomic chromosomes
  values_disomic <- colMeans(avg_chromo_matrix[disomic_chromosomes, cells_to_select])
  #calculate the ratio between the two types of chromosomes
  ratio <- log2(values_polysomic/values_disomic)
  ratio[is.na(ratio)] <- 0
  ratio[is.infinite(ratio)] <- NA
  
  #append the ratios
  cells_meta[names(ratio),]$aneu_ratio <- ratio
}

#Now it will use kmeans to determine if cells are aneuploid or normal
#NTS: for now it assumes that both normal and aneuploid cells co-exist in the data. Fix it later.
cells_meta$aneu_group <- kmeans(cells_meta$aneu_ratio, centers = 2, nstart = 25)$cluster
group_name <- cells_meta %>%
  group_by(aneu_group)%>%
  summarise(mean = mean(aneu_ratio)) %>%
  arrange(mean)
group_name$name <- c("baseline", "aneuploid")
cells_meta$aneu_group <- group_name$name[match(cells_meta$aneu_group, group_name$aneu_group)]

seuratObj@meta.data <- cells_meta[match(colnames(seuratObj), rownames(cells_meta)),]

#plotting the result
seuratObj@misc$colors$aneu_group <- c(aneuploid = "#226185", baseline = "#e6d99a")






