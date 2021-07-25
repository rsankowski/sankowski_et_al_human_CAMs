library(Seurat)
library(SingleCellExperiment)
library(clustree)
library(tidyverse)
library(scDblFinder)
library(assertthat)
library(DescTools)

#load data
load(file.path("data", "men_cp_seurat_10x.RData"))

cams <- subset(all, idents=c("9","10"))

varfeat <- 10000
cams <- SCTransform(cams, variable.features.n = varfeat)

cams <- RunPCA(cams, verbose = FALSE)
ElbowPlot(cams)
cams <- RunUMAP(cams, reduction = "pca", dims = 1:20)

#cluster data
cams <- FindNeighbors(cams, dims = 1:20)
cams <- FindClusters(cams, resolution = seq(from = .2, to = 1.6, by = .2))

clustree(cams)

#final clustering
cams <- FindClusters(cams, resolution = .8)

DimPlot(cams, label = T)
DimPlot(cams, group.by = "Compartment", label = T)

#plot gene expression
FeaturePlot(cams, "MRC1")
FeaturePlot(cams, "HLA-DRA")
FeaturePlot(cams, "P2RY12")
FeaturePlot(cams, "S100A11")
FeaturePlot(cams, "VIM")
FeaturePlot(cams, "CCR2")

#determine the most common cell type per group

save(cams, file = file.path("data", "men_cp_cams_seurat_10x.RData"))
