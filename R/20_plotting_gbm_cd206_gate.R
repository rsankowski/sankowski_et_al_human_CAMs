library(Seurat)
library(tidyverse)
library(assertthat)
library(ggpubr)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(tidyquant)
library(RColorBrewer)

source(file.path("R", "functions.R"))

load(file.path("data","GBM-SCORE-integration.RData"))

#load additional information from the indexing data
load(file.path("data","Index_anon","index_all_gbm.RData"))

#see how many complete cases are present and see that the information for only approximately 10 cells is missing
sum(complete.cases(index_all[colnames(all),])) 

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
    colSums(all[["SCT"]]@data[na.omit(signature_genes$macrophages),]) >2 ~ "BAMs",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$cDC2)[na.omit(signature_genes$cDC2) %in% rownames(all)],]) > 2 ~ "cDC2",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$cDC1),]) >2 ~ "cDC1",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$MigDCs),]) >2 ~ "MigDCs",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$tcell),]) >2 | all$celltype %in% c("CD4 TCM","CD8 TEM","Treg","CD4 TEM","gdT","CD8 TCM","CD4 TCM","CD4 Proliferating") ~ "T cells",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$pDC)[na.omit(signature_genes$pDC) %in% rownames(all)],]) > 2 ~ "pDCs",
    all$celltype %in% c("Eryth","HSPC","B intermediate","NK") ~ "other",
    T ~ all$celltype
) , levels = c("BAMs","MG","cDC2","cDC1","pDC","CD14 Mono","CD16 Mono","T cells","Oligodendr.","other")
)

#add missing values based on azimuth assignments
all$celltype_cor[is.na(all$celltype_cor)] <-  case_when(all$celltype[is.na(all$celltype_cor)] %in% c("Eryth","HSPC","B intermediate","NK") ~ "other",
            all$celltype[is.na(all$celltype_cor)] %in% c("CD4 TCM","CD8 TEM","Treg","CD4 TEM") ~ "T cells",
            T ~ all$celltype[is.na(all$celltype_cor)])

#set cluster order
order_clusters <- data.frame(seurat_clusters= all@meta.data[,"seurat_clusters"], row.names = rownames(all@meta.data)) %>%
  bind_cols(as.data.frame(t(all[["SCT"]]@scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_all(.funs=mean) %>%
  as.data.frame()

rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder the clusters
levels(all) <- order_clusters[c(7:9,5,6,10,11,4,3,1,2)]#[c(9,10,7,11,8,6,4,5,3,1,2)]

#extract metadata
metadata <- all@meta.data
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = levels(all))

#plot compartments
all %>% 
  DimPlot(group.by = "Compartment_comp", pt.size = 4) +
    scale_color_brewer(palette = "Set1") +
    theme_void()

ggsave(file.path("plots", "umap", "celseq_gbm", "compartments.pdf"), useDingbats=F)

all %>% 
  DimPlot(group.by = "Compartment_comp", pt.size = 4) +
  scale_color_brewer(palette = "Set1") +
  theme_void() +
  NoLegend()

ggsave(file.path("plots", "umap", "celseq_gbm", "compartments_no_legends.pdf"), useDingbats=F)

#marimekko compartments
all[[]][,c("seurat_clusters", "Compartment","Compartment_comp","Diagnosis")] %>% 
  na.omit() %>% 
  mutate(seurat_clusters=factor(seurat_clusters, levels = order_clusters[c(7:9,5,6,10,11,4,3,1,2)])) %>% 
mosaicGG2("seurat_clusters", "Compartment_comp") +
  scale_fill_brewer(palette = "Set1")

ggsave(file.path("plots","others","celseq_gbm","clusters_compartment_marimekko.pdf"), useDingbats=F)

all[[]][,c("seurat_clusters", "Compartment","Compartment_comp","Diagnosis")] %>% 
  na.omit() %>% 
  hyper_test_n(var1="seurat_clusters", var2="Compartment_comp") %>% 
  write_csv(file.path("data", "analysis_output_cd206_gbm", "cluster_compartment_hyper_test.csv"))

#plot celltypes 
DimPlot(all, group.by = "celltype_cor", pt.size = 4, label=T) +
  scale_color_tq() +
  theme_void()

ggsave(file.path("plots", "umap", "celseq_gbm", "celltypes.pdf"), useDingbats=F)

DimPlot(all, group.by = "celltype_cor", pt.size = 4) +
  scale_color_tq() +
  theme_void() +
  NoLegend()

ggsave(file.path("plots", "umap", "celseq_gbm", "celltypes_no_legend.pdf"), useDingbats=F)

#marimekko celltypes
na.omit(metadata) %>% 
  mosaicGG2("seurat_clusters", "celltype_cor") +
  scale_fill_tq()
ggsave(file.path("plots","others","celseq_gbm","clusters_celltype_cor.pdf"), useDingbats=F)

metadata %>% 
  na.omit() %>% 
  hyper_test_n(var1="seurat_clusters", var2="celltype_cor") %>% 
  write_csv(file.path("data", "analysis_output_cd206_gbm", "cluster_celltype_cor_hyper_test.csv"))

#plot diagnosis 
DimPlot(all, group.by = "Diagnosis", pt.size = 4, label=T) +
  scale_color_brewer(palette = "Set2") +
  theme_void()

ggsave(file.path("plots", "umap", "celseq_gbm", "celltypes.pdf"), useDingbats=F)

DimPlot(all, group.by = "Diagnosis", pt.size = 4) +
  scale_color_brewer(palette = "Set2") +
  theme_void() +
  NoLegend()

ggsave(file.path("plots", "umap", "celseq_gbm", "celltypes_no_legend.pdf"), useDingbats=F)

#marimekko celltypes
metadata[,c("seurat_clusters","Diagnosis")] %>% 
  na.omit() %>% 
  mosaicGG2("seurat_clusters", "Diagnosis") +
  scale_fill_brewer(palette = "Set2") 
  
ggsave(file.path("plots","others","celseq_gbm","clusters_Diagnosis.pdf"), useDingbats=F)

metadata[,c("seurat_clusters","Diagnosis")] %>% 
  na.omit() %>% 
  hyper_test_n(var1="seurat_clusters", var2="Diagnosis") %>% 
  write_csv(file.path("data", "analysis_output_cd206_gbm", "cluster_Diagnosis_hyper_test.csv"))

#plot clusters
DimPlot(all,pt.size = 4, label = T) +
  scale_color_manual(values = brewer.pal(3, "Pastel2")) +
  theme_void() +
  labs(title="Clusters")

ggsave(file.path("plots", "umap", "celseq_gbm", "clusters.pdf"), useDingbats=F) # images were saved in 8.68 x 5.73

DimPlot(all,pt.size = 4, label = T) +
  scale_color_manual(values = brewer.pal(3, "Pastel2")) +
  theme_void() +
  labs(title="Clusters") +
  NoLegend()

ggsave(file.path("plots", "umap", "celseq_gbm", "clusters_no_legend.pdf"), useDingbats=F) # images were saved in 8.68 x 5.73

#Heatmaps
if (!file.exists(file.path("data", "analysis_output_cd206_gbm","diffgenes_cd206_ctrl.csv"))) {
  
  all.markers<-FindAllMarkers(all,only.pos=F,min.pct=.2,loMenc.threshold=.2,return.thresh = 0.05) #,only.pos=T,min.pct=.25,loMenc.threshold=.25,return.thresh = 0.05
  
  save(all.markers, file = file.path("data","analysis_output_cd206_gbm","diffgenes_cd206_ctrl.RData"))
  write_csv(all.markers, file.path("data", "analysis_output_cd206_gbm","diffgenes_cd206_ctrl.csv"))
  
} else {
  all.markers <- read_csv( file.path("data", "analysis_output_cd206_gbm","diffgenes_cd206_ctrl.csv"))
}

all.markers <- all.markers[!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", all.markers$gene),]

top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 

heat <- DoHeatmap(all,features = top10$gene, group.colors = brewer.pal(3, "Pastel2"))
heat + scale_fill_viridis(option = "A")

ggsave(file.path("plots","heatmaps","celseq_gbm","top10-gene-heatmap-viridis-A.pdf"), width = 20, height = 16)

DoHeatmap(all,features = top10$gene, group.colors = brewer.pal(3, "Pastel2"), group.bar = F)+
  theme_void() +
  scale_fill_viridis(option = "A")  + 
  NoLegend() 

ggsave(file.path("plots","heatmaps","celseq_gbm","top10-gene-heatmap-viridis-A-heatmap-only.png"), width = 12, height = 8)


#top 25 genes
top25 <- all.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC) 

heat <- DoHeatmap(all,features = top25$gene, group.colors = brewer.pal(3, "Pastel2"))
heat + scale_fill_viridis(option = "A")

ggsave(file.path("plots","heatmaps","celseq_gbm","top25-gene-heatmap-viridis-A.pdf"), width = 20, height = 16)

DoHeatmap(all,features = top25$gene, group.colors = brewer.pal(3, "Pastel2"), group.bar = F)+
  theme_void() +
  scale_fill_viridis(option = "A")  + 
  NoLegend() 

ggsave(file.path("plots","heatmaps","celseq_gbm","top25-gene-heatmap-viridis-A-heatmap-only.png"), width = 12, height = 8)

#plot cell signatures
for (i in colnames(signature_genes)) {
  tryCatch({
    plt <- plot_expmap_seurat(na.omit(signature_genes[[i]]), point_size = 4, object = all, .retain_cl = levels(all), reduction = "NetUMAP") + labs(subtitle= paste0(i,' Signature'))
    print(plt)
    ggsave(file.path("plots", "umap", "celseq_gbm", paste0(i,"_signature.pdf")), useDingbats=F)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#run side by side analyses of the cells
all2 <- all
Idents(all2) <- paste(all$Diagnosis, all$celltype_cor, sep="_")

#volcano plots
if (!file.exists(file.path("data", "analysis_output_cd206_gbm", "all_diffgenes_Ctrl_BAMs_GBM_BAMs.RData"))) {
  all_bams_genes <- FindMarkers(all2, 
                                ident.1 = "Ctrl_BAMs",
                                ident.2 = "GBM_BAMs",
                                loMenc.threshold = 0.01,
                                min.pct = 0.01) 
  
  save(all_bams_genes, file = file.path("data", "analysis_output_cd206_gbm", "all_diffgenes_Ctrl_BAMs_GBM_BAMs.RData"))
} else {
  load(file.path("data", "analysis_output_cd206_gbm", "all_diffgenes_Ctrl_BAMs_GBM_BAMs.RData"))
}

#export suppl table 
write.csv(all_bams_genes, file.path("data", "analysis_output_cd206_gbm","Table_S_volcano_Ctrl_BAMs_GBM_BAMs.csv"))

#remove unannotated genes
all_bams_genes <- all_bams_genes %>%
  rownames_to_column(var="gene") %>% 
  filter(!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", .$gene))

top_ctrl_gbm <- all_bams_genes %>%
  filter(p_val_adj < .05) %>% 
  mutate(.rank=-log(p_val_adj)*avg_log2FC,
         genotype=case_when(avg_log2FC>0 ~ "Ctrl",
                            T~"GBM")) %>% #create a pseudo value to sort differentially expressed genes
  group_by(genotype) %>% 
  #modified from url: https://stackoverflow.com/questions/44111161/change-variable-value-for-the-first-row-group-by-subject-id-using-dplyr
  mutate(show_genes = ifelse(row_number() < 16, gene, NA)) %>% #
  dplyr::select(gene, show_genes)

all_bams_genes <- all_bams_genes %>%
  left_join(top_ctrl_gbm) %>%
  mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_log2FC) > .25, "sig.", "not sig."))


all_volcano <- ggplot(all_bams_genes, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  annotate("rect", xmin=-Inf,xmax=0,ymin=-Inf,ymax=Inf,
           fill=brewer.pal(3, "Pastel2")[2])+
  annotate("rect" ,xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf,
           fill=brewer.pal(3, "Pastel2")[1]) +
  geom_point(size=5) + 
  geom_text_repel(size=10, box.padding=1.15, max.overlaps =300, fontface="italic") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. log2FC", y="-log10 transf. adj. p-value") +
  expand_limits(x=c(-2.5,2.5))

all_volcano 

ggsave(file.path("plots", "others", "celseq_gbm", "all_Ctrl_BAMs_GBM_BAMs_volcano.pdf"), useDingbats=F, width = 11, height = 12)

#ctrl mg vs gbm mg
if (!file.exists(file.path("data", "analysis_output_cd206_gbm", "all_diffgenes_Ctrl_MG_GBM_MG.RData"))) {
  all_mg_genes <- FindMarkers(all2, 
                                ident.1 = "Ctrl_MG",
                                ident.2 = "GBM_MG",
                                loMenc.threshold = 0.01,
                                min.pct = 0.01) 
  
  save(all_mg_genes, file = file.path("data", "analysis_output_cd206_gbm", "all_diffgenes_Ctrl_MG_GBM_MG.RData"))
} else {
  load(file.path("data", "analysis_output_cd206_gbm", "all_diffgenes_Ctrl_MG_GBM_MG.RData"))
}

#export suppl table 
write.csv(all_mg_genes, file.path("data", "analysis_output_cd206_gbm","Table_S_volcano_Ctrl_MG_GBM_MG.csv"))

#remove unannotated genes
all_mg_genes <- all_mg_genes %>%
  rownames_to_column(var="gene") %>% 
  filter(!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", .$gene))

top_ctrl_gbm <- all_mg_genes %>%
  filter(p_val_adj < .05) %>% 
  mutate(.rank=-log(p_val_adj)*avg_log2FC,
         genotype=case_when(avg_log2FC>0 ~ "Ctrl",
                            T~"GBM")) %>% #create a pseudo value to sort differentially expressed genes
  group_by(genotype) %>% 
  #modified from url: https://stackoverflow.com/questions/44111161/change-variable-value-for-the-first-row-group-by-subject-id-using-dplyr
  mutate(show_genes = ifelse(row_number() < 16, gene, NA)) %>% #
  dplyr::select(gene, show_genes)

all_bams_genes <- all_mg_genes %>%
  left_join(top_ctrl_gbm) %>%
  mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_log2FC) > .25, "sig.", "not sig."))

all_volcano <- ggplot(all_bams_genes, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  annotate("rect", xmin=-Inf,xmax=0,ymin=-Inf,ymax=Inf,
           fill=brewer.pal(3, "Pastel2")[2])+
  annotate("rect" ,xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf,
           fill=brewer.pal(3, "Pastel2")[1]) +
  geom_point(size=5) + 
  geom_text_repel(size=10, box.padding=1.15, max.overlaps = 50, fontface="italic") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. log2FC", y="-log10 transf. adj. p-value") +
  expand_limits(x=c(-2.5,2.5))

all_volcano 

ggsave(file.path("plots", "others", "celseq_gbm", "all_Ctrl_MG_GBM_MG_volcano.pdf"), useDingbats=F, width = 11, height = 12)

#ctrl cDC2 vs gbm cDC2
if (!file.exists(file.path("data", "analysis_output_cd206_gbm", "all_diffgenes_Ctrl_cDC2_GBM_cDC2.RData"))) {
  all_cDC2_genes <- FindMarkers(all2, 
                                ident.1 = "Ctrl_cDC2",
                                ident.2 = "GBM_cDC2",
                                loMenc.threshold = 0.01,
                                min.pct = 0.01) 
  
  save(all_cDC2_genes, file = file.path("data", "analysis_output_cd206_gbm", "all_diffgenes_Ctrl_cDC2_GBM_cDC2.RData"))
} else {
  load(file.path("data", "analysis_output_cd206_gbm", "all_diffgenes_Ctrl_cDC2_GBM_cDC2.RData"))
}

#export suppl table 
write.csv(all_cDC2_genes, file.path("data", "analysis_output_cd206_gbm","Table_S_volcano_Ctrl_cDC2_GBM_cDC2.csv"))

#remove unannotated genes
all_cDC2_genes <- all_cDC2_genes %>%
  rownames_to_column(var="gene") %>% 
  filter(!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", .$gene))

top_ctrl_gbm <- all_cDC2_genes %>%
  filter(p_val_adj < .05) %>% 
  mutate(.rank=-log(p_val_adj)*avg_log2FC,
         genotype=case_when(avg_log2FC>0 ~ "Ctrl",
                            T~"GBM")) %>% #create a pseudo value to sort differentially expressed genes
  group_by(genotype) %>% 
  #modified from url: https://stackoverflow.com/questions/44111161/change-variable-value-for-the-first-row-group-by-subject-id-using-dplyr
  mutate(show_genes = ifelse(row_number() < 16, gene, NA)) %>% #
  dplyr::select(gene, show_genes)

all_bams_genes <- all_cDC2_genes %>%
  left_join(top_ctrl_gbm) %>%
  mutate(genes_sig = ifelse(.$p_val_adj < .05 & abs(.$avg_log2FC) > .25, "sig.", "not sig."))

all_volcano <- ggplot(all_bams_genes, aes(x=avg_log2FC, y= -log10(p_val_adj), label=show_genes, color=genes_sig)) +
  annotate("rect", xmin=-Inf,xmax=0,ymin=-Inf,ymax=Inf,
           fill=brewer.pal(3, "Pastel2")[2])+
  annotate("rect" ,xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf,
           fill=brewer.pal(3, "Pastel2")[1]) +
  geom_point(size=5) + 
  geom_text_repel(size=10, box.padding=1.15, max.overlaps = 50, fontface="italic") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        text = element_text(size=25)) +
  scale_color_manual(values = c("light grey", "black")) +
  NoLegend() +
  labs(x="avg. log2FC", y="-log10 transf. adj. p-value") +
  expand_limits(x=c(-2.5,2.5))

all_volcano 

ggsave(file.path("plots", "others", "celseq_gbm", "all_Ctrl_cDC2_GBM_cDC2_volcano.pdf"), useDingbats=F, width = 11, height = 12)




#other plots


#donut plot of cell types
meta2 <- metadata %>%
  filter(!Pat_ID %in% c("DBL-0",NA)) %>% 
  group_by(Pat_ID, Compartment, celltype_cor) %>%
  summarise(freq=n()) %>% 
  mutate(rel_freq = freq / sum(freq))

meta3 <- meta2 %>% 
  group_by(Compartment, celltype_cor) %>% 
  summarise(mean_n = mean(freq),
            std = sd(freq),
            sem = sd(freq)/sqrt(n()))


#plot donuts
meta2 %>%
  ggplot(aes(x=2, y=freq,fill=celltype_cor)) +
  geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
  coord_polar(theta='y', start=0) +
  theme_void() +
  scale_fill_manual(values=unname(glasbey)[-1]) +
  facet_wrap(~Compartment) +
  xlim(0.5, 2.5)

ggsave(file.path("plots","others","celseq_gbm","celltypes_donut.pdf"), useDingbats=F)

meta3 %>%
  ggplot(aes(x=2, y=mean_n,fill=celltype_cor)) +
  geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
  coord_polar(theta='y', start=0) +
  theme_void() +
  scale_fill_manual(values=unname(glasbey)[-1]) +
  facet_wrap(~Compartment) +
  xlim(0.5, 2.5)

ggsave(file.path("plots","others","celseq_gbm","celltypes_donut_means.pdf"), useDingbats=F)

