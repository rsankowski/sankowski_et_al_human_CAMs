library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
keytypes(org.Hs.eg.db)
library(viridis)
library(pheatmap)
library(readxl)
library(Seurat)
library(stringr)

#load data
load(file.path("data", "men_cp_cams_seurat_10x.RData"))

source(file.path("R", "functions.R"))

cams.markers <- read.csv( file.path("data","analysis_output_10x_cams","diffgenes_10x_cams.csv"))
rownames(cams.markers) <- cams.markers$gene

enrich_up <- go_term_analysis_seurat(.df=cams.markers, 
                                     .sc=cams)

write.csv(enrich_up, file.path("data","analysis_output_10x_cams","GO_terms_cams", "bp", "bp_GO_terms.csv"))

#mf terms
enrich_up_mf <- go_term_analysis_seurat(.df=cams.markers, 
                                 .sc=cams,
                                 ontogeny = 'MF')

write.csv(enrich_up_mf, file.path("data","analysis_output_10x_cams","GO_terms_cams", "mf", "mf_GO_terms.csv"))

#load filtered terms
enrich_up <- read_excel(file.path("data" ,"selcted_bp_GO_terms_for_figure.xlsx"))
terms <- unique(enrich_up$Description)

enrich_up <- read_csv(file.path("data","GO_terms_cams", "bp", "bp_GO_terms.csv")) %>%
  filter(Description %in% terms)

#set cluster order
order_clusters <- data.frame(seurat_clusters= cams@meta.data[,"seurat_clusters"], row.names = rownames(cams@meta.data)) %>%
  bind_cols(as.data.frame(t(cams[["SCT"]]@scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_cams(.funs=mean) %>%
  as.data.frame()

rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

ord_clust <- order_clusters[c(10:8,15, 7, 13:14, 16:17, 12, 11, 3:6, 1:2)]
enrich_up$Cluster <- factor(enrich_up$Cluster, levels = ord_clust)
enrich_up$Description <- reorder(enrich_up$Description,enrich_up$Description,FUN=length)

#re-order enrich_up based on the levels of Cluster
enrich_up <- enrich_up[with(enrich_up, order(Cluster)),] #from url: https://stackoverflow.com/questions/1296646/how-to-sort-a-dataframe-by-columns
enrich_up$Description <- str_to_sentence(enrich_up$Description, locale = "en")

enrich_up$Description <- factor(enrich_up$Description, levels = rev(enrich_up$Description[!duplicated(enrich_up$Description)]))
colnames(enrich_up)[10] <- 'GeneCount'

enrich_up$p.adjust <- as.numeric(enrich_up$p.adjust)

enrich_up <- enrich_up %>% 
  mutate(Celltype = factor(
    case_when(
      Cluster %in% c("7", "9", "1", "13", "4", "3", "6") ~ "CAMs",
      Cluster %in% c("0", "5", "10", "2") ~ "MG",
      Cluster == "8" ~ "Ly6clo",
      Cluster == "11" ~ "Ly6chi",
      Cluster == "12" ~ "DCs",
      Cluster == "16" ~ "Lymphocytes",
      Cluster == "15" ~ "Stromal",
      Cluster == "14" ~ "Prolif."
    ), levels = c("CAMs", "MG", "Ly6clo", "Ly6chi", "DCs", "Lymphocytes", "Stromal", "Prolif.")
  )
  )
assert_that(sum(is.na(enrich_up$Celltype))==0)

dot_plot <- ggplot(enrich_up, aes(Cluster, Description, size = GeneCount, fill= -log2(p.adjust))) + #[enrich_up$GeneCount>4,]
  geom_point(pch=21, stroke=0.25) +
  scale_fill_viridis('-log2 transf. \nadj. p-value)') +
  theme_light() +
  theme(text=element_text(size=10),
        axis.title.y=element_blank())
dot_plot

ggsave(file.path("plots", "others", "cams", 'top40-PC_cams_GoTerm_dot_plot_new.pdf'), height = 7.5, width = 7, units = 'in', useDingbats=F)

dot_plot <- ggplot(enrich_up, aes(Cluster, Description, size = GeneCount, fill= -log2(p.adjust))) + #[enrich_up$GeneCount>4,]
  geom_point(pch=21, stroke=0.25) +
  scale_fill_viridis('-log2 transf. \nadj. p-value)') +
  theme_light() +
  theme(text=element_text(size=10),
        axis.title.y=element_blank()) +
  facet_grid(. ~ Celltype,
             space = "free",
             scales = "free")
dot_plot

ggsave(file.path("plots", "others", "cams", 'top40-PC_cams_GoTerm_dot_plot_faceted.pdf'), height = 10, width = 10, units = 'in', useDingbats=F)
