library(Seurat)
library(tidyverse)
library(assertthat)
library(Polychrome)
library(ggpubr)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(ggridges)
library(fishualize)

source(file.path("R", "functions.R"))

load(file.path("data","SCORE-integration.RData"))
data("glasbey")

#load additional information from the indexing data
load(file.path("data","Index_anon","index_all_ctrl.RData")) 

#see how many complete cases are present and see that the information for only approximately 10 cells is missing
sum(complete.cases(index_all[colnames(all),])) 

#assign metadata
all$Compartment <- index_all[colnames(all),"Compartment_ind"]
all$Pat_ID <- index_all[colnames(all),"Pat_ID"]

#assign celltype information based on gene expression and clustering results
all$celltype_cor <- factor(
  case_when(
    colSums(all[["SCT"]]@data[na.omit(signature_genes$monocytes),]) >2 ~ "CD14 Mono",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$microglia),]) >2 & all$seurat_clusters == "4" ~ "MG",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$macrophages),]) >2 ~ "BAMs",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$cDC2)[na.omit(signature_genes$cDC2) %in% rownames(all)],]) >2 ~ "cDC2",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$cDC1),]) >2 ~ "cDC1",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$MigDCs),]) >2 ~ "MigDCs",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$oligodendrocyte),]) >2 ~ "Oligodendr.",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$tcell),]) >2 | all$celltype %in% c("CD4 TCM","CD8 TEM","Treg","CD4 TEM") ~ "T cells",
    colSums(all[["SCT"]]@data[na.omit(signature_genes$pDC)[na.omit(signature_genes$DC) %in% rownames(all)],]) >2 ~ "pDCs",
    all$celltype %in% c("Eryth","HSPC","B intermediate","NK") ~ "other",
    T ~ all$celltype
) , levels = c("BAMs","MG","cDC2","cDC1","pDC","CD14 Mono","CD16 Mono","T cells","other")
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
levels(all) <- order_clusters[c(9,10,7,11,8,6,4,5,3,1,2)]

#extract metadata
metadata <- all@meta.data
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = levels(all))

#plot compartments
all %>% 
na.omit() %>% 
  droplevels() %>% 
  DimPlot(group.by = "Compartment", pt.size = 4) +
    scale_color_brewer(palette = "Set1") +
    theme_void()

ggsave(file.path("plots", "umap", "celseq_cd206", "compartments.pdf"), useDingbats=F)

all %>% 
  na.omit() %>% 
  droplevels() %>% 
DimPlot(group.by = "Compartment", pt.size = 4) +
  scale_color_brewer(palette = "Set1") +
  theme_void() +
  NoLegend()

ggsave(file.path("plots", "umap", "celseq_cd206", "compartments_no_legends.pdf"), useDingbats=F)

#marimekko compartments
na.omit(metadata) %>% 
mosaicGG2("seurat_clusters", "Compartment") +
  scale_fill_brewer(palette = "Set1")

ggsave(file.path("plots","others","celseq_cd206","clusters_compartment_marimekko.pdf"), useDingbats=F)

metadata %>% 
  na.omit() %>% 
  hyper_test_n(var1="seurat_clusters", var2="Compartment") %>% 
    write_csv(file.path("data", "analysis_output_10x", "cluster_compartment_hyper_test.csv"))

#plot celltypes
#define colors 
cols <- c(colorRampPalette(brewer.pal(8, "Set2"))(8),colors_fig[3])

DimPlot(all, group.by = "celltype_cor", pt.size = 4, label=T) +
  #scale_color_manual(values = unname(glasbey)[-1]) +
  #scale_color_fish(option = "Centropyge_loricula", discrete = T) +
  scale_color_manual(values = cols) +
  theme_void()

ggsave(file.path("plots", "umap", "celseq_cd206", "celltypes.pdf"), useDingbats=F)

DimPlot(all, group.by = "celltype_cor", pt.size = 4) +
  scale_color_manual(values = cols) +
  theme_void() +
  NoLegend()

ggsave(file.path("plots", "umap", "celseq_cd206", "celltypes_no_legend.pdf"), useDingbats=F)

#marimekko celltypes
na.omit(metadata) %>% 
  mosaicGG2("seurat_clusters", "celltype_cor") +
  #scale_fill_manual(values = colors_fig2)
scale_fill_manual(values = cols)
ggsave(file.path("plots","others","celseq_cd206","clusters_celltype_cor.pdf"), useDingbats=F)

metadata %>% 
  na.omit() %>% 
  hyper_test_n(var1="seurat_clusters", var2="celltype_cor") %>% 
  write_csv(file.path("data", "analysis_output_10x", "cluster_celltype_cor_hyper_test.csv"))

#plot clusters
DimPlot(all,pt.size = 4, label = T) +
  scale_color_manual(values = c(colors_pat, colors_many, colors_fig)[-2]) +
  theme_void() +
  labs(title="Clusters")

ggsave(file.path("plots", "umap", "celseq_cd206", "clusters.pdf"), useDingbats=F) # images were saved in 8.68 x 5.73

DimPlot(all,pt.size = 4, label = T) +
  scale_color_manual(values = c(colors_pat, colors_many, colors_fig)[-2]) +
  theme_void() +
  labs(title="Clusters") +
  NoLegend()

ggsave(file.path("plots", "umap", "celseq_cd206", "clusters_no_legend.pdf"), useDingbats=F) # images were saved in 8.68 x 5.73

#Heatmaps
if (!file.exists(file.path("data", "analysis_output_CD206_gate","diffgenes_cd206_ctrl.csv"))) {
  
  all.markers<-FindAllMarkers(all,only.pos=F,min.pct=.2,loMenc.threshold=.2,return.thresh = 0.05) #,only.pos=T,min.pct=.25,loMenc.threshold=.25,return.thresh = 0.05
  
  save(all.markers, file = file.path("data","analysis_output_CD206_gate","diffgenes_cd206_ctrl.RData"))
  write_csv(all.markers, file.path("data", "analysis_output_CD206_gate","diffgenes_cd206_ctrl.csv"))
  
} else {
  all.markers <- read_csv( file.path("data", "analysis_output_CD206_gate","diffgenes_cd206_ctrl.csv"))
}

all.markers <- all.markers[!grepl("^(HTRA|LIN|EEF|CTC-|MIR|CTD-|AC0|RP|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS|PMAIP1|HSP|NEAT1|HIST|MALAT1|RP)", all.markers$gene),]

top10 <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 

heat <- DoHeatmap(all,features = top10$gene, group.colors = c(colors_pat, colors_many, colors_fig)[-2])
heat + scale_fill_viridis(option = "A")

ggsave(file.path("plots","heatmaps","celseq_cd206","top10-gene-heatmap-viridis-A.pdf"), width = 20, height = 16)

DoHeatmap(all,features = top10$gene, group.colors = c(colors_pat, colors_many, colors_fig)[-2], group.bar = F)+
  theme_void() +
  scale_fill_viridis(option = "A")  + 
  NoLegend() 

ggsave(file.path("plots","heatmaps","celseq_cd206","top10-gene-heatmap-viridis-A-heatmap-only.png"), width = 12, height = 8)


#top 25 genes
top25 <- all.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC) 

heat <- DoHeatmap(all,features = top25$gene, group.colors = c(colors_pat, colors_many, colors_fig)[-2])
heat + scale_fill_viridis(option = "A")

ggsave(file.path("plots","heatmaps","celseq_cd206","top25-gene-heatmap-viridis-A.pdf"), width = 20, height = 16)

DoHeatmap(all,features = top25$gene, group.colors = c(colors_pat, colors_many, colors_fig)[-2], group.bar = F)+
  theme_void() +
  scale_fill_viridis(option = "A")  + 
  NoLegend() 

ggsave(file.path("plots","heatmaps","celseq_cd206","top25-gene-heatmap-viridis-A-heatmap-only.png"), width = 12, height = 8)

#plot cell signatures
for (i in colnames(signature_genes)) {
  tryCatch({
    plt <- plot_expmap_seurat(na.omit(signature_genes[[i]]), point_size = 4, object = all, .retain_cl = levels(all), reduction = "NetUMAP") + labs(subtitle= paste0(i,' Signature'))
    print(plt)
    ggsave(file.path("plots", "umap", "celseq_cd206", paste0(i,"_signature.pdf")), useDingbats=F)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

#plot individual genes
walk(c("SIGLEC1","LYVE1","HLA-DRA","MRC1","CD163","STAB1"), function(x) {
  plot_expmap_seurat(x, point_size = 4, object = all, .retain_cl = levels(all), reduction = "NetUMAP")
  ggsave(file.path("plots", "umap", "celseq_cd206", paste0(x,".pdf")), useDingbats=F)
})

all2 <- subset(all, subset = celltype_cor=="BAMs")
Idents(all2) <- all2$Compartment

#prepare didge plots of the data
genes <- c("SIGLEC1","LYVE1","HLA-DRA","MRC1","CD163")

all3 <- all2[["SCT"]]@data[genes,] %>% 
  t %>% 
  as.matrix %>% 
  as.data.frame %>% 
  bind_cols(all2[[]][,c("Compartment","seurat_clusters")]) %>% 
  na.omit %>% 
  rownames_to_column(var = "cell_ID") %>% 
  pivot_longer(SIGLEC1:CD163, names_to = "gene", values_to = "abundance") %>% 
  mutate(gene=factor(.$gene, levels=c("MRC1","LYVE1","CD163","SIGLEC1","HLA-DRA")))

#ridge plots
all3 %>% 
  ggplot(aes(x = abundance, y = Compartment, fill = Compartment)) +
    geom_density_ridges() +
    scale_fill_brewer(palette = "Set1") +
    facet_wrap(~gene, ncol=1) +
    theme_pubclean() + 
    theme(legend.position = "none",
          text=element_text(size=60)) 

ggsave(file.path("plots", "others", "celseq_cd206", "single_gene_ridgeplot.pdf"), useDingbats=F, height = 25.95, width = 11.01)
  
#compare index and gene expression data
genes2 <- c("MRC1","HLA-DRA","CX3CR1","ITGAM","CD14","PTPRC")
all_gene <- all2[["SCT"]]@data[genes2,] %>% 
  t %>% 
  as.matrix %>% 
  as.data.frame %>% 
  bind_cols(all2[[]][,c("Compartment","seurat_clusters")]) %>% 
  na.omit %>% 
  rownames_to_column(var = "cell_ID") %>% 
  pivot_longer(MRC1:PTPRC, names_to = "gene", values_to = "abundance") %>% 
  group_by(gene) %>% 
  mutate(norm_abundance=(abundance-(min(abundance)))/(max(abundance)-min(abundance)),
         molecule="Transcript",
         gene=str_replace_all(gene,pattern = c("MRC1"="MRC1/CD206","ITGAM"="ITGAM/CD11b","PTPRC"="PTPRC/CD45")),
         gene=factor(gene, levels=c("MRC1/CD206","HLA-DRA","CX3CR1","ITGAM/CD11b","CD14","PTPRC/CD45")))

prots <- index_all[unique(all3$cell_ID),c("cell_ID","CD206","MHCII","Cx3Cr1","CD11b","CD14","CD45")] %>% 
  pivot_longer(CD206:CD45, names_to="gene", values_to="abundance") %>% 
  group_by(gene, .drop=F) %>% 
  mutate(norm_abundance=(abundance-min(abundance))/(max(abundance)-min(abundance)),
         molecule="Protein",
         gene=str_replace_all(gene,pattern = c("CD206"="MRC1/CD206","MHCII"="HLA-DRA","Cx3Cr1"="CX3CR1","CD11b"="ITGAM/CD11b","CD45"="PTPRC/CD45")),
         gene=factor(gene, levels=c("MRC1/CD206","HLA-DRA","CX3CR1","ITGAM/CD11b","CD14","PTPRC/CD45")))

#check if both datasets have equal numbers of rows
assert_that(nrow(all_gene)==nrow(prots))

#combine both datasets
both <- bind_rows(all_gene, prots)

#calculate the number of 0 values in the transcript data
pct_zero <- both %>% 
  group_by(molecule, gene) %>% 
  summarise(percent_zero=sum(abundance==0)/length(abundance) * 100)

both %>% ggplot(aes(x=norm_abundance,color=molecule)) +
  stat_ecdf(size=1) +
  facet_wrap(~gene, scales = "free") +
  geom_text(data=pct_zero[pct_zero$molecule=="Transcript",], aes(x=0.0, y=.95,label=paste0(round(percent_zero,1)," %")), size=5, nudge_x=.1,inherit.aes = F) +
  theme_pubclean() +
  theme(text=element_text(size=15),
        legend.position = "bottom") +
  scale_y_continuous(name = "Proportion of cells",breaks = c(0,.5,1), labels = c("0","0.5","1")) +
  scale_x_continuous(name = "Min-Max-normalized molecule abundance",breaks = c(0,.5,1), labels = c("0","0.5","1")) +
  scale_color_fish_d(option = "Cirrhilabrus_solorensis",direction = 1,begin = .1, end=.9) +
  expand_limits(x=c(0,1))
  
ggsave(file.path("plots", "others", "celseq_cd206", "ecdf_plots_transcriptome_index_data.pdf"), useDingbats=F)

#plot violin plots of genes and cell cycle state
markers <- c("S.Score","G2M.Score","KCNQ1OT1")

walk(markers, function(x){
  plt <- VlnPlot(all, x) +
    scale_fill_manual(values = c(colors_pat, colors_many, colors_fig)[-2]) +
    theme_pubclean() +
    theme(text=element_text(size=30),
          legend.position = "none") +
    labs(x="Cluster")
  print(plt)
  
  ggsave(file.path("plots", "umap", "celseq_cd206", paste0(x,"_violin_plot.pdf")), useDingbats=F, height = 11, width = 11)
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

ggsave(file.path("plots", "others", "celseq_cd206", "all_c9_c10_volcano.pdf"), useDingbats=F, width = 11, height = 12)

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

ggsave(file.path("plots","others","celseq_cd206","celltypes_donut.pdf"), useDingbats=F)

meta3 %>%
  ggplot(aes(x=2, y=mean_n,fill=celltype_cor)) +
  geom_bar(position = 'fill', stat = 'identity', color='black', lwd=0.1) +
  coord_polar(theta='y', start=0) +
  theme_void() +
  scale_fill_manual(values=unname(glasbey)[-1]) +
  facet_wrap(~Compartment) +
  xlim(0.5, 2.5)

ggsave(file.path("plots","others","celseq_cd206","celltypes_donut_means.pdf"), useDingbats=F)

