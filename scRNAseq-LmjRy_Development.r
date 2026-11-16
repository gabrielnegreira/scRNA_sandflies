library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(gdata)
library(rmarkdown)
library(limma)
library(SingleR)
library(sctransform)
library(devtools)
library(forcats)
library(DT)
library(data.table)
library(RColorBrewer)

set.seed(42)

# scRNAseq analysis - Development dataset

#Load datasets----

#LmjRy-AMA
Ama <- Read10X(data.dir = '~/bioinfo/LmjRy_scRNA_Development_Transm/Amastigotes-Aggr/filtered_feature_bc_matrix')
#LmjRy-Pro_Day2 infection
Day2 <- Read10X(data.dir = '~/bioinfo/LmjRy_scRNA_Development_Transm/Day2-Aggr-infection/filtered_feature_bc_matrix')
#LmjRy-Pro_Day5 infection
Day5 <- Read10X(data.dir = '~/bioinfo/LmjRy_scRNA_Development_Transm/Day5-Aggr-infection/filtered_feature_bc_matrix')
#LmjRy-Pro_Day9 infection
Day9 <- Read10X(data.dir = '~/bioinfo/LmjRy_scRNA_Development_Transm/Day9-Aggr-infection/filtered_feature_bc_matrix')

#Create Seurat objects----

#LmjRy-Ama_Day0 infection
LmjRyAMA <- CreateSeuratObject(counts = Ama, min.cells = 2, min.features = 200, project = "AMA")
LmjRyDay2 <- CreateSeuratObject(counts = Day2, min.cells = 2, min.features = 200, project = "Day2-Inf")
LmjRyDay5 <- CreateSeuratObject(counts = Day5, min.cells = 2, min.features = 200, project = "Day5-Inf")
LmjRyDay9 <- CreateSeuratObject(counts = Day9, min.cells = 2, min.features = 200, project = "Day9-Inf")

#Regressing mitochondrial genes----

#LmjRy-Ama_Day0 infection
LmjRyAMA[["percent.mt"]]<-PercentageFeatureSet(LmjRyAMA,pattern="LmjF.maxi")
VlnPlot(LmjRyAMA, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3, col="lightblue")
LmjRyAMA<-subset(LmjRyAMA, subset= nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<20)
VlnPlot(LmjRyAMA, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3, col="lightblue")

#LmjRy-Pro_Day2 infection
LmjRyDay2[["percent.mt"]]<-PercentageFeatureSet(LmjRyDay2,pattern="LmjF.maxi")
VlnPlot(LmjRyDay2, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3, col="lightblue")
LmjRyDay2<-subset(LmjRyDay2, subset= nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
VlnPlot(LmjRyDay2, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3, col="lightblue")

#LmjRy-Pro_Day5 infection
LmjRyDay5[["percent.mt"]]<-PercentageFeatureSet(LmjRyDay5,pattern="LmjF.maxi")
VlnPlot(LmjRyDay5, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3, col="lightblue")
LmjRyDay5<-subset(LmjRyDay5, subset= nFeature_RNA>100 & nFeature_RNA<5000 & percent.mt<10)
VlnPlot(LmjRyDay5, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3, col="lightblue")

#LmjRy-Pro_Day9 infection
LmjRyDay9[["percent.mt"]]<-PercentageFeatureSet(LmjRyDay9,pattern="LmjF.maxi")
VlnPlot(LmjRyDay9, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3, col="lightblue")
LmjRyDay9<-subset(LmjRyDay9, subset= nFeature_RNA>100 & nFeature_RNA<5000 & percent.mt<20)
VlnPlot(LmjRyDay9, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3, col="lightblue")

#SCTranform each sample separately for SCT assay----

Ama.sct <- SCTransform(LmjRyAMA, vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = T)
Day2.sct <- SCTransform(LmjRyDay2, vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = T)
Day5.sct <- SCTransform(LmjRyDay5, vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = T)
Day9.sct <- SCTransform(LmjRyDay9, vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = T)

#Normalize data for RNA assay----

DefaultAssay(Ama.sct) <- "RNA"
DefaultAssay(Day2.sct) <- "RNA"
DefaultAssay(Day5.sct) <- "RNA"
DefaultAssay(Day9.sct) <- "RNA"

#Amastigote
Ama.sct <-NormalizeData(Ama.sct,normalization.method="LogNormalize",scale.factor=10000)
Ama.sct <-FindVariableFeatures(Ama.sct, selection.method = "vst", nfeatures = 2000)
#Day2
Day2.sct <-NormalizeData(Day2.sct,normalization.method="LogNormalize",scale.factor=10000)
Day2.sct <-FindVariableFeatures(Day2.sct, selection.method = "vst", nfeatures = 2000)
#Day5
Day5.sct <-NormalizeData(Day5.sct,normalization.method="LogNormalize",scale.factor=10000)
Day5.sct <-FindVariableFeatures(Day5.sct, selection.method = "vst", nfeatures = 2000)
#Day9
Day9.sct <-NormalizeData(Day9.sct,normalization.method="LogNormalize",scale.factor=10000)
Day9.sct <-FindVariableFeatures(Day9.sct, selection.method = "vst", nfeatures = 2000)

#Dimentional reduction for each sample----

#Amastigotes
DefaultAssay(Ama.sct) <- "SCT"
Ama.sct <- ScaleData(Ama.sct)
Ama.sct <- RunPCA(Ama.sct, npcs = 30)
ElbowPlot(Ama.sct, ndims = 30, reduction ="pca") 
Ama.sct <- RunUMAP(Ama.sct, dims = 1:14)
Ama.sct <- RunTSNE(Ama.sct, dims = 1:14)
Ama.sct <- FindNeighbors(Ama.sct, dims = 1:14)

#Day2
DefaultAssay(Day2.sct) <- "SCT"
Day2.sct <- ScaleData(Day2.sct)
Day2.sct <- RunPCA(Day2.sct, npcs = 30)
ElbowPlot(Day2.sct, ndims = 30, reduction ="pca") 
Day2.sct <- RunUMAP(Day2.sct, dims = 1:14)
Day2.sct <- RunTSNE(Day2.sct, dims = 1:14)
Day2.sct <- FindNeighbors(Day2.sct, dims = 1:14)


#Day5
DefaultAssay(Day5.sct) <- "SCT"
Day5.sct <- ScaleData(Day5.sct)
Day5.sct <- RunPCA(Day5.sct, npcs = 30)
ElbowPlot(Day5.sct, ndims = 30, reduction ="pca") 
Day5.sct <- RunUMAP(Day5.sct, dims = 1:14)
Day5.sct <- RunTSNE(Day5.sct, dims = 1:14)
Day5.sct <- FindNeighbors(Day5.sct, dims = 1:14)


#Day9
DefaultAssay(Day9.sct) <- "SCT"
Day9.sct <- ScaleData(Day9.sct)
Day9.sct <- RunPCA(Day9.sct, npcs = 30)
ElbowPlot(Day9.sct, ndims = 30, reduction ="pca") 
Day9.sct <- RunUMAP(Day9.sct, dims = 1:14)
Day9.sct <- RunTSNE(Day9.sct, dims = 1:14)
Day9.sct <- FindNeighbors(Day9.sct, dims = 1:14)

#Merge Seurat objects----

LmjRyDevelopment <- merge(Ama.sct, y = c(Day2.sct, Day5.sct, Day9.sct), add.cell.ids = c("Ama", "Day2", "Day5", "Day9"), project = "LmjRyDevelopment",
                      merge.data = TRUE)

#Plot UMAP----

VariableFeatures(LmjRyDevelopment[["SCT"]]) <- rownames(LmjRyDevelopment[["SCT"]]@scale.data) 

LmjRyDevelopment <- RunPCA(LmjRyDevelopment, npcs = 30)
ElbowPlot(LmjRyDevelopment, ndims = 30, reduction ="pca") 
LmjRyDevelopment <- RunUMAP(LmjRyDevelopment, dims = 1:14)
LmjRyDevelopment <- RunTSNE(LmjRyDevelopment, dims = 1:14)
LmjRyDevelopment <- FindNeighbors(LmjRyDevelopment, dims = 1:14)
LmjRyDevelopment <- FindClusters(LmjRyDevelopment, resolution = 0.14) #6 clusters

#Find subclusters from cluster 4 and 5----

DefaultAssay(LmjRyDevelopment) <- "SCT"

LmjRyDevelopment <- FindSubCluster(LmjRyDevelopment, cluster = "2", graph.name = "SCT_snn", subcluster.name = "sub2.clusters", resolution = 0.1)
Idents(LmjRyDevelopment) <- LmjRyDevelopment$sub2.clusters
LmjRyDevelopment <- FindSubCluster(LmjRyDevelopment, cluster = "4", graph.name = "SCT_snn", subcluster.name = "sub5.sub4.clusters", resolution = 0.05)
Idents(LmjRyDevelopment) <- LmjRyDevelopment$sub5.sub4.clusters

LmjRyDevelopment <- RenameIdents(LmjRyDevelopment, `1` = "1", `2` = "2", `3` = "3",
                             `4_1` = "4_1", `4_0` = "4_2", `5_1` = "5_1", `5_0` = "5_2", `6` = "6")