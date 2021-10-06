library(Seurat)
library(clustree)
library(tidyverse)
library(assertthat)
library(Polychrome)
library(fishualize)
library(ggpubr)
library(ggrepel)
library(tidyquant)

source(file.path("R", "functions.R"))

load(file.path("data", "men_cp_seurat_10x.RData"))
data("glasbey")

#set cluster order
order_clusters <- data.frame(seurat_clusters= all@meta.data[,"seurat_clusters"], row.names = rownames(all@meta.data)) %>%
  bind_cols(as.data.frame(t(all[["SCT"]]@scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_all(.funs=mean) %>%
  as.data.frame()

rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder the clusters
levels(all) <- order_clusters[c(11,15,16,18,17,13,14,1,5:10,2:4,12)]

#extract metadata
metadata <- all@meta.data
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = levels(all))

#plot compartments
DimPlot(all, group.by = "Compartment", pt.size = 2) +
  scale_color_brewer(palette = "Set1") +
  theme_void()

ggsave(file.path("plots", "umap", "10x", "compartments.pdf"), useDingbats=F)

DimPlot(all, group.by = "Compartment", pt.size = 2) +
  scale_color_brewer(palette = "Set1") +
  theme_void() +
  NoLegend()

ggsave(file.path("plots", "umap", "10x", "compartments_no_legends.pdf"), useDingbats=F)

#marimekko compartments
mosaicGG2(metadata, "seurat_clusters", "Compartment") +
  scale_fill_brewer(palette = "Set1")

ggsave(file.path("plots","others","10x","clusters_compartment_marimekko.pdf"), useDingbats=F)

hyper_test_n(metadata, var1="seurat_clusters", var2="Compartment") %>% 
  write_csv(file.path("data", "analysis_output_10x", "cluster_compartment_hyper_test.csv"))

#plot celltypes
DimPlot(all, group.by = "Celltype", pt.size = 2, label=T) +
  #scale_color_manual(values = unname(glasbey)[-1]) +
  #scale_color_fish(option = "Centropyge_loricula", discrete = T) +
  scale_color_manual(values = c(colors_many[-c(21,22)], colors_fig)) +
  theme_void()

ggsave(file.path("plots", "umap", "10x", "celltypes.pdf"), useDingbats=F)

DimPlot(all, group.by = "Celltype", pt.size = 2) +
  scale_color_manual(values = unname(glasbey)[-1]) +
  theme_void() +
  NoLegend()

ggsave(file.path("plots", "umap", "10x", "celltypes_no_legend.pdf"), useDingbats=F)

#plot clusters
DimPlot(all,pt.size = 2, label = T) +
  scale_color_manual(values = c(colors_pat, colors_many, colors_fig)[-2]) +
  theme_void() +
  labs(title="Clusters")

ggsave(file.path("plots", "umap", "10x", "clusters.pdf"), useDingbats=F) # images were saved in 8.68 x 5.73

DimPlot(all,pt.size = 2, label = T) +
  scale_color_manual(values = c(colors_pat, colors_many, colors_fig)[-2]) +
  theme_void() +
  labs(title="Clusters") +
  NoLegend()

ggsave(file.path("plots", "umap", "10x", "clusters_no_legend.pdf"), useDingbats=F) # images were saved in 8.68 x 5.73

#donut plot of cell types
meta2 <- metadata %>%
  filter(!Patient_ID %in% c("DBL-0",NA)) %>% 
  group_by(Patient_ID, Compartment, Cell_lineage,Celltype) %>%
  summarise(freq=n()) %>% 
  mutate(rel_freq = freq / sum(freq))

meta3 <- meta2 %>% 
  group_by(Compartment, Cell_lineage,Celltype) %>% 
  summarise(mean_n = mean(freq),
            std = sd(freq),
            sem = sd(freq)/sqrt(n()))
  

#plot donuts
meta2 %>%
  ggplot(aes(x=2, y=freq,fill=Celltype)) +
  geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
  coord_polar(theta='y', start=0) +
  theme_void() +
  scale_fill_manual(values=unname(glasbey)[-1]) +
  facet_wrap(~Compartment) +
  xlim(0.5, 2.5)

ggsave(file.path("plots","others","10x","celltypes_donut.pdf"), useDingbats=F)

meta3 %>%
  group_by(Compartment, Cell_lineage) %>% 
  summarise(cum_n = sum(mean_n)) %>% 
  ggplot(aes(x=2, y=cum_n,fill=Cell_lineage)) +
  geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
  coord_polar(theta='y', start=0) +
  theme_void() +
  scale_fill_tq() +
  facet_wrap(~Compartment) +
  xlim(0.5, 2.5)

ggsave(file.path("plots","others","10x","celltypes_donut_means.pdf"), useDingbats=F)

meta4 %>% 
  ggplot(aes(Compartment, mean_n, fill=Compartment)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0) +
  geom_point(pch=21, size=5) +
  facet_wrap(~Cell_lineage, nrow=1) +
  stat_compare_means(method = "t.test")

meta2 %>% 
  ggplot(aes(Compartment, freq, fill=Compartment)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0) +
  geom_point(pch=21, size=5) +
  facet_wrap(~Celltype, scales = "free") +
  stat_compare_means()

#stat testing
meta4 <- meta2 %>% 
  group_by(Patient_ID, Compartment, Cell_lineage) %>% 
  summarise(mean_n = mean(freq))#
wilcox.test(mean_n ~ Compartment + Cell_lineage, data = meta4)


meta3 %>%
  group_by(Compartment, Cell_lineage) %>% 
  summarise(cum_n = sum(mean_n)) %>% 
  hyper_test_n(var1 = "Compartment",var2 = "Cell_lineage")

#marimekko plots celltypes
mosaicGG2(metadata, "Celltype", "Compartment") +
  scale_fill_brewer(palette = "Set1")

ggsave(file.path("plots","others","10x","celltype_compartment_marimekko.pdf"), useDingbats=F)

hyper_test_n(metadata, var1="Celltype", var2="Compartment") %>% 
  write_csv(file.path("data", "analysis_output_10x", "celltype_compartment_hyper_test.csv"))

#plot celltypes separated by cell lineages
walk(unique(metadata$Cell_lineage), function(x){
  .df <- droplevels(metadata[metadata$Cell_lineage == x, ]) %>% 
    left_join(dplyr::count(.,Celltype)) %>% 
    mutate(Celltype=reorder(.$Celltype, -n))
    
  plt <- mosaicGG2(.df, "Celltype", "Compartment") +
    scale_fill_brewer(palette = "Set1")
  print(plt)
  
  ggsave(file.path("plots","others","10x",paste(x,"faceted_celltype_compartment_marimekko.pdf", sep="_")), useDingbats=F)
  
})

#Heatmaps
if (!file.exists( file.path("data", "analysis_output_10x","diffgenes_10x.csv"))) {
  
  all.markers<-FindAllMarkers(all,only.pos=F,min.pct=.2,loMenc.threshold=.2,return.thresh = 0.05) #,only.pos=T,min.pct=.25,loMenc.threshold=.25,return.thresh = 0.05
  
  save(all.markers, file = file.path("data", "analysis_output_10x","diffgenes_10x.RData"))
  write_csv(all.markers, file.path("data", "analysis_output_10x","diffgenes_10x.csv"))
  
} else {
  all.markers <- read_csv( file.path("data", "analysis_output_10x","diffgenes_10x.csv"))
}

all.markers <- all.markers[!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", all.markers$gene),]

top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 

heat <- DoHeatmap(all,features = top10$gene, group.colors = c(colors_pat, colors_many, colors_fig)[-2])
heat + scale_fill_viridis(option = "A")

ggsave(file.path("plots","heatmaps","10x","top10-gene-heatmap-viridis-A.pdf"), width = 20, height = 16)

DoHeatmap(all,features = top10$gene, group.colors = c(colors_pat, colors_many, colors_fig)[-2], group.bar = F)+
  theme_void() +
  scale_fill_viridis(option = "A")  + 
  NoLegend() 

ggsave(file.path("plots","heatmaps","10x","top10-gene-heatmap-viridis-A-heatmap-only.png"), width = 12, height = 8)

all.markers <- all.markers[!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", all.markers$gene),]

top20 <- all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 

heat <- DoHeatmap(all,features = top20$gene, group.colors = c(colors_pat, colors_many, colors_fig)[-2])
heat + scale_fill_viridis(option = "A")

ggsave(file.path("plots","heatmaps","10x","top20-gene-heatmap-viridis-A.pdf"), width = 20, height = 16)

DoHeatmap(all,features = top20$gene, group.colors = c(colors_pat, colors_many, colors_fig)[-2], group.bar = F)+
  theme_void() +
  scale_fill_viridis(option = "A")  + 
  NoLegend() 

ggsave(file.path("plots","heatmaps","10x","top20-gene-heatmap-viridis-A-heatmap-only.png"), width = 12, height = 8)

#plot cell signatures
signature_genes <- data.frame("monocytes"=c('CCR2', 'CLEC12A', 'PLAC8', 'FCN1', 'S100A9'),
                              "macrophages"=c("MRC1", "MS4A7", "CD163", "LYVE1", "STAB1"),
                              "microglia"= c('P2RY12', 'CX3CR1', 'CSF1R', 'TMEM119', 'SLC2A5'),
                              "tcell"=c('TRAC', 'TRBC2', 'CD52', 'IL32', NA),
                              "myeloid"=c('ITGAM',  'MS4A6A', 'TYROBP', 'CD14', NA),
                              "oligodendrocyte"=c('MBP',  'MOG', 'MAG', 'PLP1', NA),
                              "bcells"=c('CD79A', 'IGHG4', 'IGLL5', NA, NA),
                              "astrocyte"=c("MenAP", "HEPACAM","SOX9","AQP4",NA),
                              "apc"=c("CD74", "CD80", "CD86", "HLA-DRA", "CD40"),
                              "pDC"=c("IL3RA","LILRA4","TCF4","SELL","LTB"),
                              "cDC1"=c("XCR1", "CLEC9A","CADM1", "IRF8","BATF3"),
                              "cDC2"=c("FCER1A", "CLEC10A", "CD1C","CST7","CCR6"),
                              "MigDCs"= c("CCR7", "LAMP3", "SAMSN1",NA,NA),
                              "Endothelial Cells"=c("HSPG2","PLVAP","FLT1","VWF","CD34"),
                              stringsAsFactors = F)

for (i in colnames(signature_genes)) {
  plt <- plot_expmap_seurat(na.omit(signature_genes[[i]]), point_size = 2, object = all, .retain_cl = levels(all)) + labs(subtitle= paste0(i,' Signature'))
  print(plt)
  ggsave(file.path("plots", "umap", "10x", paste0(i,"_signature_cp.pdf")), useDingbats=F)
}

#plot s score and g2m score
walk(c("S.Score","G2M.Score"), function(x) {
  FeaturePlot(all, x, pt.size = 2) +
  scale_color_gradientn(colors = c("darkblue","lightblue2","yellow","red2")) +
  theme_void()
  
  ggsave(file.path("plots", "umap", "10x", paste0(x,".pdf")), useDingbats=F)
  })

#volcano plots
#compare clusters 9 and 10
if (!file.exists(file.path("data", "analysis_output_10x", "all_diffgenes_cl9_cl10.RData"))) {
  all_c9c10genes <- FindMarkers(all, 
                                  ident.1 = "9",
                                  ident.2 = "10",
                                  loMenc.threshold = 0.01,
                                  min.pct = 0.01) 
  
  save(all_c9c10genes, file = file.path("data", "analysis_output_10x", "all_diffgenes_cl9_cl10.RData"))
} else {
  load(file.path("data", "analysis_output_10x", "all_diffgenes_cl9_cl10.RData"))
}

#export suppl table 
write.csv(all_c9c10genes, file.path("data", "analysis_output_10x","Table_S_volcano_10x_cl9_cl10.csv"))

#remove unannotated genes
all_c9c10genes <- all_c9c10genes %>%
  rownames_to_column(var="gene") %>% 
  filter(!grepl("(Gm|Rp|Rik|mt-|RP)", .$gene))

top5_c9 <- all_c9c10genes %>%
  mutate(.rank=-log(p_val_adj)*avg_log2FC) %>% #create a pseudo value to sort differentially expressed genes
  top_n(15, .rank) %>%
  mutate(show_genes = gene) %>%
  dplyr::select(gene, show_genes)

top_5_both <- all_c9c10genes %>%
  mutate(.rank=-log(p_val_adj)*avg_log2FC) %>% #create a pseudo value to sort differentially expressed genes
  top_n(-15, .rank) %>%
  mutate(show_genes = gene) %>%
  dplyr::select(gene, show_genes) %>%
  bind_rows(top5_c9) 

all_c9c10genes <- all_c9c10genes %>%
  left_join(top_5_both) %>%
  mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_log2FC) > .25, "sig.", "not sig."))


all_volcano <- ggplot(all_c9c10genes, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  annotate("rect", xmin=-Inf,xmax=0,ymin=-Inf,ymax=Inf,
           fill=c(colors_pat, colors_many, colors_fig)[-2][2], alpha =.4)+
  annotate("rect" ,xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf,
           fill=c(colors_pat, colors_many, colors_fig)[-2][1], alpha =.4) +
  geom_point(size=5) + 
  geom_text_repel(size=10, box.padding=1.15, max.overlaps = 20, fontface="italic") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. log2FC", y="-log10 transf. adj. p-value") +
  expand_limits(y=35, x=c(-12,12))

all_volcano 

ggsave(file.path("plots", "others", "10x", "all_c9_c10_volcano.pdf"), useDingbats=F, width = 11, height = 12)



#plot violin plots of individual genes
#create an overview of the expression of selected genes
all2 <- all %>% 
  subset(idents= c("9","10"))

#make sure that the cell ids are identical
assert_that(identical(rownames(all2[[]]), names(all2[["SCT"]]@data["P2RY12",])))

single_genes <- bind_cols(all2[[]], as.data.frame(t(as.matrix(all2[["SCT"]]@data[c("P2RY12","S100A11","VIM"),])))) %>%
  pivot_longer(P2RY12:VIM, names_to = "gene", values_to = "expression") 

single_genes %>% 
  ggplot(aes(gene, expression,fill=Compartment)) +
  geom_violin(scale = "width") +
  theme_pubclean()











#other plots

#marimekko plots cell lineages
mosaicGG2(metadata, "Cell_lineage", "Compartment") +
  scale_fill_brewer(palette = "Set1")

ggsave(file.path("plots","others","10x","cell_lineage_compartment_marimekko.pdf"), useDingbats=F)

hyper_test_n(metadata, var1="Cell_lineage", var2="Compartment") %>% 
  write_csv(file.path("data", "analysis_output_10x", "cell_lineage_compartment_hyper_test.csv"))

#plot dot plots
walk(levels(meta2$Cell_lineage), function(x) {
  a <- meta2 %>% 
    filter(Cell_lineage==x) %>% 
    ggplot(aes(Compartment, freq, fill=Celltype, group=Celltype)) +
    geom_point(size=5, pch=21) +
    stat_summary(fun = mean, geom="crossbar") +
    geom_line(data = meta3[meta3$Cell_lineage==x, ], aes(Compartment, y=mean_n)) +
    scale_fill_manual(values=unname(glasbey)[-1]) +
    #facet_grid(Cell_lineage~Celltype,scales="free",space="free") +
    facet_wrap(~Celltype, nrow=1)+
    expand_limits(y=0) +
    theme_pubclean() 
  
  print(a)
  
  ggsave(file.path("plots","others","10x",paste("celltype_dotplot",x,".pdf",sep="_")), useDingbats=F, height = 6, width=12)
  
})

#bray-curtis dissimilarity of leptomeningeal and choroid plexus macrophages
set.seed(79106)

#copare based on variable genes
a_cp <- all[["SCT"]]@counts[VariableFeatures(all),] %>%
  as.matrix() %>%
  as.data.frame() %>%
  dplyr::select(colnames(all)[all$seurat_clusters %in% c("9","10") & all$Compartment == "CP"]) %>%
  t() 

a_cp <- a_cp %>%
  vegdist(method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
          na.rm = FALSE) %>%
  as.numeric()

a_men <- all[["SCT"]]@counts[VariableFeatures(all),] %>%
  as.matrix() %>%
  as.data.frame() %>%
  dplyr::select(colnames(all)[all$seurat_clusters %in% c("9","10") & all$Compartment == "Men"]) %>%
  t() 

a_men <- a_men %>%
  vegdist(method="bray", binary=FALSE, diag=FALSE, upper=FALSE,
          na.rm = FALSE) %>%
  as.numeric()

both <- bind_rows(data.frame(condition="Men", dist=a_men),
                 data.frame(condition="CP", dist=a_cp))
both$condition <- factor(both$condition, levels = c("CP", "Men"))

both %>%
  ggplot(aes(condition, dist, fill=condition)) +
  geom_violin(draw_quantiles = 0.5) +
  scale_fill_manual(values = colors_many[2:4]) +
  stat_compare_means() +
  theme_linedraw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=20)) +
  labs(title = "MG Bray-Curtis Dissimilarity", y="Bray-Curtis Dissimilarity Coefficient")

ggsave(file.path("plots", "others", "all", "all_bray_curtis_violin_plots.pdf"), useDingbats=F)

#comparison between leptomeningeal and choroid plexus macrophages

all2 <- all
Idents(all2) <- paste(all2$Compartment,all2$Celltype, sep = "_")

if (!file.exists(file.path("data", "analysis_output_10x", "all_diffgenes_cp_vs_Men.RData"))) {
  all_men_cp_genes <- FindMarkers(all2, 
                                  ident.1 = "Men_CAMs",
                                  ident.2 = "CP_CAMs",
                                  loMenc.threshold = 0.01,
                                  min.pct = 0.01)
  
  save(all_men_cp_genes, file = file.path("data", "analysis_output_10x", "all_diffgenes_cp_vs_Men.RData"))
} else {
  load(file.path("data", "analysis_output_10x", "all_diffgenes_cp_vs_Men.RData"))
}

#export suppl table 8
write.csv(all_men_cp_genes, file.path("data", "analysis_output_10x", "Table_all_men_cp_diffgenes.csv"))
#remove unannotated genes
all_men_cp_genes <- all_men_cp_genes %>%
  rownames_to_column(var="gene") %>% 
  filter(!grepl("(Gm|Rp|Rik|mt-|RP)", .$gene))

top15_men <- all_men_cp_genes %>%
  mutate(.rank=-log(p_val_adj)*avg_log2FC) %>% #create a pseudo value to sort differentially expressed genes
  top_n(15, .rank) %>%
  mutate(show_genes = gene) %>%
  dplyr::select(gene, show_genes)

top_15_both <- all_men_cp_genes %>%
  mutate(.rank=-log(p_val_adj)*avg_log2FC) %>% #create a pseudo value to sort differentially expressed genes
  top_n(-15, .rank) %>%
  mutate(show_genes = gene) %>%
  dplyr::select(gene, show_genes) %>%
  bind_rows(top15_men) 

all_men_cp_genes <- all_men_cp_genes %>%
  left_join(top_15_both) %>%
  mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_log2FC) > .25, "sig.", "not sig."))

all_men_cp_volcano <- ggplot(all_men_cp_genes, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  annotate("rect", xmin=-Inf,xmax=0,ymin=-Inf,ymax=Inf,
           fill=colors_many[2], alpha =.4)+
  annotate("rect" ,xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf,
           fill=colors_many[3], alpha =.4) +
  geom_point(size=5) + 
  geom_text_repel(size=7, box.padding=1.15, max.overlaps = 20) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. log2FC", y="-log10 transf. adj. p-value") 

all_men_cp_volcano 

ggsave(file.path("plots", "others", "all", "men_cp_volcano.pdf"), useDingbats=F)

#export data
write.csv(tibble(Cell_ID=rownames(all[[]]), UMAP_1 = all@reductions$umap@cell.embeddings[,1], UMAP_2 = all@reductions$umap@cell.embeddings[,2], all[[]][,c(2:7,18:24)]), file.path("data", "analysis_output_10x", "Table_S1_metadata_10x.csv"), row.names = F)
