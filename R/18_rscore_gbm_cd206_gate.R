#packages were installed according to https://github.com/wycwycpku/RSCORE

library(tidyverse)
library(Seurat)
library(Matrix)
library(cowplot)
library(clustree)
library(RSCORE)
library(igraph)
library(genesorteR)
library(harmony)
library(cellity)

date <- Sys.Date()

#load data
load(file.path("data","prdata_celseq_adult_cd206_gate_only.RData"))

  #create seurat object and add metadata
  all <- CreateSeuratObject(prdata, min.cells = 5, min.features = 500) %>% 
    AddMetaData(object=.,
                metadata= gsub("_.*", "", colnames(.)),
                col.name="ID") %>%
    #assess mitochondrial genes
    AddMetaData(object=.,
                metadata=PercentageFeatureSet(., pattern = "^MT-"),
                col.name="percent.mt") %>%
    #assess ribosomal genes
    AddMetaData(object=.,
                metadata=PercentageFeatureSet(., pattern = "^(RPS|RPL)"),
                col.name="percent.ribo") %>%#add azimuth cell type assignment run on 2021-07-08 on https://app.azimuth.hubmapconsortium.org/app/human-pbmc
    AddMetaData(object=.,
                metadata=read.delim(file.path("data", "azimuth_cd206_pbmc", "azimuth_pred.tsv"), row.names = "cell")[colnames(.),"predicted.celltype.l2"],
                col.name="celltype") %>% 
    AddMetaData(object=.,
                metadata=ifelse(grepl("GBM", colnames(.)),"GBM", "Control"),
                col.name="Diagnosis") %>% 
    #remove cells with more than 25% mt genes and more than 5000 genes
    subset(subset = percent.mt <= 25 & nFeature_RNA <= 5000) 
  
  #preprocess seurat object
  all <- SCTransform(all, vars.to.regress = "percent.mt", verbose = FALSE,variable.features.n=10000)
 
  #run RSCORE
  hs_network <- as.matrix(readRDS(system.file('extdata','9606_ppi_matrix_BioGRID-3.5.173.Rda',package = 'RSCORE')))
  all <- R.SCORE(Data = all, PPI = hs_network, max_step = 10, nCores = 4)
  all <- RunPCA(all, features = rownames(all), assay = "Net", npcs = 50,
                       reduction.name = "NetPCA", reduction.key = "NetPCA_", verbose = F)
  
  ElbowPlot(all, reduction = "NetPCA")
  
  all <- RunUMAP(all, reduction = "NetPCA", dims = 1:30,
                        reduction.name = "NetUMAP",  reduction.key = "NetUMAP_")
  
  DimPlot(all, reduction = 'NetUMAP', pt.size = 1, group.by = 'orig.ident')
  
  all <- FindNeighbors(all, reduction = "NetPCA",dims = 1:30)
  all <- FindClusters(all, resolution = seq(.2, 2, by=0.2))
  
  #url https://cran.r-project.org/web/packages/clustree/vignettes/clustre e.html#seurat-objects
  clustree(all)
  
  all <- FindClusters(all, resolution = 1)
  
  DimPlot(all, reduction = "NetUMAP", group.by = "celltype", label = T) #, split.by = 'batch'
  DimPlot(all, reduction = "NetUMAP", group.by = "Diagnosis", label = T) #, split.by = 'batch'
  DimPlot(all, label = T)
  
  DefaultAssay(all) <- "SCT"
  FeaturePlot(all, "HLA-DRA")
  FeaturePlot(all, "MRC1")
  FeaturePlot(all, "TMEM119")
  FeaturePlot(all, "P2RY12")
  FeaturePlot(all, "VIM")
  
  #cell cycle scoring
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  # Create our Seurat object and complete the initalization steps
  all <- CellCycleScoring(all, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  Idents(all) <- all$seurat_clusters
  
  #save harmony integration
  save(all, file = file.path("data","GBM-SCORE-integration.RData"))
