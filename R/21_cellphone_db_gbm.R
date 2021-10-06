library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(assertthat)

source(file.path("R", "functions.R"))

load(file.path("data","GBM-SCORE-integration.RData"))

#load additional information from the indexing data
load(file.path("data","Index_anon","index_all_gbm.RData"))

#assign metadata
all$Compartment <- factor(index_all[colnames(all),"Compartment_ind"], levels = c("CP","Dura","MenM","PVM"))
all$Diagnosis <- factor(index_all[colnames(all),"Diagnosis_ind"], levels = c("Ctrl","GBM"))
all$Pat_ID <- index_all[colnames(all),"Pat_ID"]
all$Compartment_comp <- factor(paste(all$Diagnosis, all$Compartment, sep="_"),levels = c("Ctrl_CP","Ctrl_Dura","Ctrl_MenM","Ctrl_PVM","GBM_PVM","GBM_MenM"))

#assign celltype information based on gene expression and clustering results
all$celltype_cor <- factor(
  case_when(
    all$seurat_clusters == "10" ~ "Oligodendr.",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$monocytes),]) >2 ~ "CD14 Mono",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$microglia),]) >2 | all$seurat_clusters == "2" ~ "MG",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$macrophages),]) >2 ~ "CAMs",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$cDC2)[na.omit(signature_genes$cDC2) %in% rownames(all)],]) > 2 ~ "cDC2",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$cDC1),]) >2 ~ "cDC1",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$MigDCs),]) >2 ~ "MigDCs",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$tcell),]) >2 | all$celltype %in% c("CD4 TCM","CD8 TEM","Treg","CD4 TEM","gdT","CD8 TCM","CD4 TCM","CD4 Proliferating") ~ "T cells",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$pDC)[na.omit(signature_genes$pDC) %in% rownames(all)],]) > 2 ~ "pDCs",
    all$celltype %in% c("Eryth","HSPC","B intermediate","NK") ~ "other",
    T ~ all$celltype
  ) , levels = c("CAMs","MG","cDC2","cDC1","pDC","CD14 Mono","CD16 Mono","T cells","Oligodendr.","other")
)

all$celltype_cor2 <- factor(
  case_when(
    all$seurat_clusters == "10" ~ "Oligodendr.",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$monocytes),]) >2 ~ "CD14 Mono",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$microglia),]) >2 | all$seurat_clusters == "2" ~ "MG",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$macrophages),]) >2 ~ "CAMs",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$cDC2)[na.omit(signature_genes$cDC2) %in% rownames(all)],]) > 2 ~ "cDC2",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$cDC1),]) >2 ~ "cDC1",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$MigDCs),]) >2 ~ "MigDCs",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$pDC)[na.omit(signature_genes$pDC) %in% rownames(all)],]) > 2 ~ "pDCs",
    all$celltype %in% c("Eryth","HSPC","B intermediate") ~ "other",
    T ~ all$celltype
  ) , levels = c("CAMs","MG","cDC2","cDC1","pDC","CD14 Mono","CD16 Mono","CD4 TCM","CD8 TEM","Treg","CD4 TEM","gdT","CD8 TCM","CD4 Proliferating","NK","Oligodendr.","other")
)
#add missing values based on azimuth assignments
all$celltype_cor2[is.na(all$celltype_cor2)] <-  case_when(all$celltype[is.na(all$celltype_cor2)] %in% c("Eryth","HSPC","B intermediate") ~ "other",
                                                        T ~ all$celltype[is.na(all$celltype_cor2)])

#make sure no values are missing
assert_that(sum(is.na(all$celltype_cor2))==0)

#extract metadata
metadata <- all@meta.data

#prepare data
counts <- all[["RNA"]]@counts
ids <- clusterProfiler::bitr(rownames(counts), fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
counts <- counts[ids$SYMBOL[!duplicated(ids$ENSEMBL)],]
rownames(counts) <- ids$ENSEMBL[!duplicated(ids$ENSEMBL)]

#export data
write.csv(counts, file.path("data","analysis_output_cd206_gbm","cellphonedb","cell_counts.csv"))
write.csv(data.frame(cell=rownames(metadata),cell_type=paste(metadata$Compartment_comp, metadata$celltype_cor2, sep = "_")), file.path("data","analysis_output_cd206_gbm","cellphonedb","cell_meta.csv"), quote = F, row.names = F)

#meningeal cells
#export data
metadata_filt <- metadata[!is.na(metadata$Compartment),]
metadata_men <- metadata_filt[which(metadata_filt$Compartment == "MenM"), ]
counts_men <- counts[,colnames(counts) %in% rownames(metadata_men)]
write.csv(counts_men, file.path("data","analysis_output_cd206_gbm","cellphonedb_men","cell_counts.csv"))
write.csv(data.frame(cell=rownames(metadata_men),cell_type=paste(metadata_men$Compartment_comp, metadata_men$celltype_cor2, sep = "_")), file.path("data","analysis_output_cd206_gbm","cellphonedb_men","cell_meta.csv"), quote = F, row.names = F)

#pvm
metadata_pvm <- metadata_filt[which(metadata_filt$Compartment == "PVM"), ]
counts_pvm <- counts[,colnames(counts) %in% rownames(metadata_pvm)]
write.csv(counts_pvm, file.path("data","analysis_output_cd206_gbm","cellphonedb_pvm","cell_counts.csv"))
write.csv(data.frame(cell=rownames(metadata_pvm),cell_type=paste(metadata_pvm$Compartment_comp, metadata_pvm$celltype_cor2, sep = "_")), file.path("data","analysis_output_cd206_gbm","cellphonedb_pvm","cell_meta.csv"), quote = F, row.names = F)

#the algorithm is run in the command line and the output is saved in the cellphone_db/out folder
#running cellphonedb using the following code from url:https://github.com/Teichlab/cellphonedb :
#for setup:
#conda create -n cpdb python=3.7
#python -m venv cpdb
#source activate cpdb
#pip install cellphonedb
#cd into the cellphonedb directory
#cellphonedb method statistical_analysis cell_meta.csv cell_counts.csv --threads=8
#cellphonedb plot heatmap_plot cell_meta.csv         
#cellphonedb plot dot_plot
#cellphonedb plot dot_plot --rows in/rows.txt --columns in/columns.txt


