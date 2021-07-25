library(Seurat)
library(clustree)
library(tidyverse)
library(assertthat)
library(Polychrome)
library(fishualize)
library(ggpubr)
library(ggrepel)
     
source(file.path("R", "functions.R"))

load(file.path("data", "men_cp_cams_seurat_10x.RData"))

#set cluster order
order_clusters <- data.frame(seurat_clusters= cams@meta.data[,"seurat_clusters"], row.names = rownames(cams@meta.data)) %>%
  bind_cols(as.data.frame(t(cams[["SCT"]]@scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_all(.funs=mean) %>%
  as.data.frame()

rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder the clusters
levels(cams) <- order_clusters

#extract metadata
metadata <- cams@meta.data
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = levels(cams))

#plot compartments
DimPlot(cams, group.by = "Compartment", pt.size = 8) +
  scale_color_brewer(palette = "Set1") +
  theme_void()

ggsave(file.path("plots", "umap", "10x_cams", "compartments.pdf"), useDingbats=F)

DimPlot(cams, group.by = "Compartment", pt.size = 8) +
  scale_color_brewer(palette = "Set1") +
  theme_void() +
  NoLegend()

ggsave(file.path("plots", "umap", "10x_cams", "compartments_no_legends.pdf"), useDingbats=F)

#marimekko compartments
mosaicGG2(metadata, "seurat_clusters", "Compartment") +
  scale_fill_brewer(palette = "Set1")

ggsave(file.path("plots","others","10x_cams","clusters_compartment_marimekko.pdf"), useDingbats=F)

hyper_test_n(metadata, var1="seurat_clusters", var2="Compartment") %>% 
  write_csv(file.path("data", "analysis_output_10x_cams", "cluster_compartment_hyper_test.csv"))

#plot clusters
DimPlot(cams,pt.size = 8, label = T) +
  scale_color_manual(values = c(colors_pat, colors_many, colors_fig)[-2]) +
  theme_void() +
  labs(title="Clusters")

ggsave(file.path("plots", "umap", "10x_cams", "clusters.pdf"), useDingbats=F) # images were saved in 8.68 x 5.73

DimPlot(cams,pt.size = 8, label = T) +
  scale_color_manual(values = c(colors_pat, colors_many, colors_fig)[-2]) +
  theme_void() +
  labs(title="Clusters") +
  NoLegend()

ggsave(file.path("plots", "umap", "10x_cams", "clusters_no_legend.pdf"), useDingbats=F) # images were saved in 8.68 x 5.73

#Heatmaps
if (!file.exists(file.path("data","analysis_output_10x_cams","diffgenes_10x_cams.csv"))) {
  
  cams.markers<-FindAllMarkers(cams,only.pos=F,min.pct=.2,loMenc.threshold=.2,return.thresh = 0.05) #,only.pos=T,min.pct=.25,loMenc.threshold=.25,return.thresh = 0.05
  
  save(cams.markers, file = file.path("data","analysis_output_10x_cams","diffgenes_10x_cams.RData"))
  write_csv(cams.markers, file.path("data","analysis_output_10x_cams","diffgenes_10x_cams.csv"))
  
} else {
  cams.markers <- read_csv( file.path("data","analysis_output_10x_cams","diffgenes_10x_cams.csv"))
}

cams.markers <- cams.markers[!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", cams.markers$gene),]

top10 <- cams.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 

heat <- DoHeatmap(cams,features = top10$gene, group.colors = c(colors_pat, colors_many, colors_fig)[-2])
heat + scale_fill_viridis(option = "A")

ggsave(file.path("plots","heatmaps","10x_cams","top10-gene-heatmap-viridis-A.pdf"), width = 20, height = 16)

DoHeatmap(cams,features = top10$gene, group.colors = c(colors_pat, colors_many, colors_fig)[-2], group.bar = F)+
  theme_void() +
  scale_fill_viridis(option = "A")  + 
  NoLegend() 

ggsave(file.path("plots","heatmaps","10x_cams","top10-gene-heatmap-viridis-A-heatmap-only.png"), width = 12, height = 8)

#plot cell signatures
signature_genes <- data.frame("monocytes"=c('CCR2', 'CLEC12A', 'PLAC8', 'FCN1', 'S100A9'),
                              "macrophages"=c("MRC1", "MS4A7", "CD163", "LYVE1", "STAB1"),
                              "microglia"= c('P2RY12', 'CX3CR1', 'CSF1R', 'TMEM119', 'SLC2A5'),
                              "tcell"=c('TRAC', 'TRBC2', 'CD52', 'IL32', NA),
                              "myeloid"=c('ITGAM',  'MS4A6A', 'TYROBP', 'CD14', NA),
                              "oligodendrocyte"=c('MBP',  'MOG', 'MAG', 'PLP1', NA),
                              "bcells"=c('CD79A', 'IGHG4', 'IGLL5', NA, NA),
                              "astrocyte"=c("MenAP", "HEPACAM","SOX9","AQP4",NA),
                              "apc"=c("CD74", "CD80", "CD86", "HLA-DRA", "CD40"), stringsAsFactors = F)

for (i in colnames(signature_genes)) {
  plt <- plot_expmap_seurat(na.omit(signature_genes[[i]]), point_size = 8, object = cams, .retain_cl = levels(cams)) + labs(subtitle= paste0(i,' Signature'))
  print(plt)
  ggsave(file.path("plots", "umap", "10x_cams", paste0(i,"_signature_cp.pdf")), useDingbats=F)
}

genes <- c("LYVE1","HLA-DRB1","MRC1","S100A11","CXCL8","FN1", "SIGLEC1", "CD163")

for (i in genes) {
  plt <- plot_expmap_seurat(i, point_size = 8, object = cams, .retain_cl = levels(cams))
  print(plt)
  ggsave(file.path("plots", "umap", "10x_cams", paste0(i,".pdf")), useDingbats=F)
}
