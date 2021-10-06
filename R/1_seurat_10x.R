library(Seurat)
library(SingleCellExperiment)
library(clustree)
library(tidyverse)
library(scDblFinder)
library(assertthat)
library(DescTools)

#load data
fls <- list.files(file.path("data", "counts_10x"))
lst <- list()
lst <- map(fls, function(i) {
  Read10X(data.dir=file.path("data", "counts_10x", i, "outs", "filtered_feature_bc_matrix")) %>% 
    CreateSeuratObject(min.cells = 5, min.features = 500) %>% 
    AddMetaData(object=.,
                metadata=PercentageFeatureSet(., pattern = "^MT-"),
                col.name="percent.mt") 
})

map(lst, dim)
walk(lst, function(i) {
  print(VlnPlot(i, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
})

#remove cells with less than 25% mitochondrial genes
lst <- map(lst, function(i) {
  subset(i, percent.mt < 25)
})

#check dataset dimensions
map(lst, dim)

names(lst) <- fls

#conduct sctransform with 10000 features
varfeat <- 10000
lst <- lapply(X = lst, FUN = SCTransform, variable.features.n = varfeat)

#sanity check of the data
assertthat::assert_that(dim(lst[[1]][["SCT"]]@scale.data)[1] == varfeat)
assertthat::assert_that(dim(lst[[2]][["SCT"]]@scale.data)[1] == varfeat)

features <- SelectIntegrationFeatures(object.list = lst, nfeatures = 10000)
lst <- PrepSCTIntegration(object.list = lst, anchor.features = features)

#integrate data
immune.anchors <- FindIntegrationAnchors(object.list = lst, normalization.method = "SCT", 
                                         anchor.features = features)
all <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

all <- RunPCA(all, verbose = FALSE)
ElbowPlot(all)
all <- RunUMAP(all, reduction = "pca", dims = 1:20)

#assign compartment
all$Compartment <- case_when(
  grepl("_1$", colnames(all)) ~ "CP",
  TRUE ~ "Men"
)

#celltypes based on azimuth human brain cells
celltypes_cp <- read_tsv(file.path("data", "azimuth_10x_brain", "cp", "azimuth_pred.tsv")) %>%
  as.data.frame()
celltypes_men <- read_tsv(file.path("data", "azimuth_10x_brain", "men", "azimuth_pred.tsv")) %>%
  as.data.frame()
rownames(celltypes_cp) <-  paste(celltypes_cp$cell, "1", sep = "_")
rownames(celltypes_men) <-  paste(celltypes_men$cell, "2", sep = "_")

all$Celltype_br <- bind_rows(celltypes_cp, celltypes_men)[colnames(all),]$predicted.subclass

#celltype assignment based on azimuth pbmc
celltypes_cp_pbmc <- read_tsv(file.path("data", "azimuth_10x_pbmc", "cp", "azimuth_pred.tsv")) %>%
  as.data.frame()
celltypes_men_pbmc <- read_tsv(file.path("data", "azimuth_10x_pbmc", "men", "azimuth_pred.tsv")) %>%
  as.data.frame()
rownames(celltypes_cp_pbmc) <-  paste(celltypes_cp_pbmc$cell, "1", sep = "_")
rownames(celltypes_men_pbmc) <-  paste(celltypes_men_pbmc$cell, "2", sep = "_")

all$Celltype_pbmc <- bind_rows(celltypes_cp_pbmc, celltypes_men_pbmc)[colnames(all),]$predicted.celltype.l2

DimPlot(all, group.by = "Compartment")
DimPlot(all, group.by = "Celltype_pbmc", label = T)
DimPlot(all, group.by = "Celltype_br")

#cluster data
all <- FindNeighbors(all, dims = 1:20)
all <- FindClusters(all, resolution = seq(from = .2, to = 1.6, by = .2))

clustree(all)

#final clustering
all <- FindClusters(all, resolution = .8)

DimPlot(all, label = T)

#plot gene expression
FeaturePlot(all, "MRC1")
FeaturePlot(all, "HLA-DRA")

#determine the most common cell type per group
most_common <- all[[]][,c("seurat_clusters","Celltype_pbmc")] %>% group_by(seurat_clusters) %>% summarise(most_common_cell_type = Mode(Celltype_pbmc)) #count() %>% group_by(seurat_clusters, Celltype_pbmc) %>% dplyr::arrange(Celltype_pbmc,desc(n)) %>% top_n(1)

#reassign CAM clusters 9,10
all$Celltype <- factor(case_when(
  all$seurat_clusters %in% c("9", "10") ~ "CAMs",
  all$Celltype_pbmc %in% c("Platelet", "Eryth") & all$seurat_clusters == "4" ~ most_common[most_common$seurat_clusters==4,]$most_common_cell_type,
  all$Celltype_pbmc %in% c("Platelet", "Eryth") & all$seurat_clusters == "5" ~ most_common[most_common$seurat_clusters==5,]$most_common_cell_type,
  all$Celltype_pbmc %in% c("Platelet", "Eryth") & all$seurat_clusters == "6" ~ most_common[most_common$seurat_clusters==6,]$most_common_cell_type,
  all$Celltype_pbmc %in% c("Platelet", "Eryth") & all$seurat_clusters == "11" ~ most_common[most_common$seurat_clusters==11,]$most_common_cell_type,
  T ~ all$Celltype_pbmc
), levels = c("CAMs","CD14 Mono","CD16 Mono","cDC1","cDC2","pDC",
            "CD4 CTL","CD4 Naive","CD4 TCM","CD4 TEM","CD8 Naive","CD8 Proliferating","CD8 TCM","CD8 TEM","dnT","gdT","ILC","MAIT","Treg",
            "B intermediate","B memory","B naive","Plasmablast",
            "NK","NK Proliferating","NK_CD56bright"))


DimPlot(all, group.by = "Celltype", label = T)

#assign celltype group
all$Cell_lineage <- factor(
  case_when(
  all$Celltype %in% c("B intermediate","B memory","B naive","Plasmablast") ~ "B cells",
  all$Celltype %in%  c("CAMs","CD14 Mono","CD16 Mono","cDC1","cDC2","pDC") ~ "Myeloid cells",
  all$Celltype %in% c("CD4 CTL","CD4 Naive","CD4 TCM","CD4 TEM","CD8 Naive","CD8 Proliferating","CD8 TCM","CD8 TEM","dnT","gdT","ILC","MAIT","Treg") ~ "T cells",
  all$Celltype %in% c("NK","NK Proliferating","NK_CD56bright") ~ "NK"
),
levels = c("Myeloid cells","T cells","B cells","NK") )

assert_that(sum(is.na(all$Cell_lineage))==0)

## Remove doublets
sce <- SingleCellExperiment(assays=list(counts=all@assays$RNA@counts))
sce <- scDblFinder(sce)
all <- all[,sce@colData$scDblFinder.class == "singlet"]

#assign patient id
cp <- read.delim(file.path("data","scSplit_result_cp.csv"),row.names=1,sep="\t") %>% 
  dplyr::rename(Patient_ID=Cluster) %>% 
  mutate(Patient_ID=gsub("SNG", "cp",.$Patient_ID))
rownames(cp) <- paste(rownames(cp),"1",sep="_")

men <- read.delim(file.path("data","scSplit_result_men.csv"),row.names = 1,sep = "\t") %>% 
  dplyr::rename(Patient_ID=Cluster) %>% 
  mutate(Patient_ID=gsub("SNG", "men",.$Patient_ID))
rownames(men) <- paste(rownames(men),"2",sep="_")

both <- rbind(cp, men)
all$Patient_ID <- both[rownames(all[[]]),]

#cell cycle scoring
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Create our Seurat object and complete the initalization steps
all <- CellCycleScoring(all, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Idents(all) <- all$seurat_clusters

save(all, file = file.path("data", "men_cp_seurat_10x.RData"))
