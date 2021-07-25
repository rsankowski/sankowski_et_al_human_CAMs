#setup the folders
dir.create("data")
dir.create(file.path("data", "counts_CD206_gate"))
dir.create(file.path("data", "counts_CD45_gate"))
dir.create(file.path("data", "counts_nn"))
dir.create(file.path("data", "counts_10x"))
dir.create(file.path("data", "counts_fetal"))
dir.create(file.path("data", "analysis_output_CD206_gate"))
dir.create(file.path("data", "analysis_output_CD45_gate"))
dir.create(file.path("data", "analysis_output_nn"))
dir.create(file.path("data", "analysis_output_10x"))
dir.create(file.path("data", "analysis_output_10x_cams"))
dir.create(file.path("data", "analysis_output_fetal"))
dir.create(file.path("data", "azimuth_10x_brain"))
dir.create(file.path("data", "azimuth_10x_brain", "cp"))
dir.create(file.path("data", "azimuth_10x_brain","men"))
dir.create(file.path("data", "azimuth_10x_pbmc"))
dir.create(file.path("data", "azimuth_10x_pbmc", "cp"))
dir.create(file.path("data", "azimuth_10x_pbmc","men"))
dir.create(file.path("data", "azimuth_cd206_brain"))
dir.create(file.path("data", "azimuth_cd45_brain"))
dir.create(file.path("data", "azimuth_cd206_pbmc"))
dir.create(file.path("data", "azimuth_cd45_pbmc"))
dir.create(file.path("data", "azimuth_gbm"))
dir.create(file.path("data", "azimuth_fetal"))
dir.create(file.path("data","analysis_output_10x_cams","GO_terms_cams"))
dir.create(file.path("data","analysis_output_10x_cams","GO_terms_cams", "bp"))
dir.create(file.path("data","analysis_output_10x_cams","GO_terms_cams", "mf"))
dir.create("R")
dir.create("plots")
dir.create(file.path("plots", "QC"))
dir.create(file.path("plots", "umap"))
dir.create(file.path("plots", "umap","10x"))
dir.create(file.path("plots", "umap","10x_cams"))
dir.create(file.path("plots", "umap","celseq_cd206"))
dir.create(file.path("plots", "umap","celseq_fetal"))
dir.create(file.path("plots", "umap","celseq_gbm"))
dir.create(file.path("plots", "heatmaps"))
dir.create(file.path("plots", "heatmaps","10x"))
dir.create(file.path("plots", "heatmaps","10x_cams"))
dir.create(file.path("plots", "heatmaps","celseq_cd206"))
dir.create(file.path("plots", "heatmaps","celseq_fetal"))
dir.create(file.path("plots", "heatmaps","celseq_gbm"))
dir.create(file.path("plots", "others"))
dir.create(file.path("plots", "others","10x"))
dir.create(file.path("plots", "others","10x_cams"))
dir.create(file.path("plots", "others","celseq_cd206"))
dir.create(file.path("plots", "others","celseq_fetal"))
dir.create(file.path("plots", "others","celseq_gbm"))

# the data object SC_NT2.rda is in the Microglia_Fillatreau folder

#install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

#the following code was adjusted from url https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
list.of.packages <- c(
  "SingleCellExperiment",
  "SingleR",
  "Seurat",
  "RColorBrewer",
  "tidyverse",
  "assertthat",
  "viridis",
  "readxl",
  "scDblFinder",
  "clustree",
  "scuttle",
  "scran",
  "AUCell", 
  "STRINGdb", 
  "propr", 
  "coop", 
  "network",
  "intergraph",
  "doRNG",
  "doParallel",
  "mahmoudibrahim/genesorteR",
  "clusterProfiler", 
  "org.Hs.eg.db",
  "scDblFinder",
  "SingleCellExperiment",
  "ggpubr",
  "DescTools",
  "ggrepel",
  "cellity"
  
)


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

#additional steps for RSCORE
install.packages("doMC", repos="http://R-Forge.R-project.org")


