library(tidyverse)
library(tidyquant)
library(Seurat)
library(clustree)
library(fishualize)
library(RColorBrewer)
library(ggpubr)

#load data
load(file.path("data","ctrl_cytof", "ctrl_cytof.RData"))
source(file.path("R","functions.R"))

#assign celltype
ctrl_cytof$Celltype <- ctrl_cytof$ClusterID
levels(ctrl_cytof$Celltype) <- c("Microglia", "other", "Microglia", "Lymphocytes", "NK cells", "Lymphocytes", "CAMs", "Lymphocytes", "Monocytes", "Lymphocytes", "DCs")

ctrl_cytof$Region <- factor(gsub("PC/PV", "PV/PC", ctrl_cytof$Region), levels = c("CP","Dura","Men","PV/PC"))
  
df2 <- ctrl_cytof %>%
group_by(ClusterID) %>%
  summarise(median_x = median(X),
            median_y = median(Y))

#tidy up protein names
colnames(ctrl_cytof)[c(8:43)] <- toupper(colnames(ctrl_cytof[,c(8:43)]))
colnames(ctrl_cytof) <- gsub("P2Y12", "P2RY12",colnames(ctrl_cytof))

#order clusters based on similarity
ord <- ctrl_cytof[,c(6,8:43)] %>% 
  na.omit() %>% 
  group_by(ClusterID) %>% 
  summarise_all(.funs = mean) %>% 
  as.data.frame()

rownames(ord) <- ord$ClusterID
ord$ClusterID <- NULL

ord_2 <- ord %>% 
  as.matrix() %>% 
  scale() %>% 
  dist() %>% 
  hclust

ctrl_cytof$ClusterID <- factor(ctrl_cytof$ClusterID, levels=rev(ord_2$order))

#Cluster
    ctrl_cytof %>%
      ggplot(aes(X, Y, color = ClusterID)) +
      geom_point(size=2) +
      theme_void() +
      geom_text(data = df2, aes(median_x, median_y, label = ClusterID), color="black", size = 10) +
      scale_color_tq() +
      guides(colour = guide_legend(override.aes = list(size=5))) #from url: https://stackoverflow.com/questions/20415963/how-to-increase-the-size-of-points-in-legend-of-ggplot2
    
    ggsave(file.path("plots","umap","ctrl_cytof","clusters_with_legend.pdf"), width = 7.7, height = 5.17, useDingbats = F)
    
    ctrl_cytof %>%
      ggplot(aes(X, Y, color = ClusterID)) +
      geom_point(size=2) +
      theme_void() +
      geom_text(data = df2, aes(median_x, median_y, label = ClusterID), color="black", size = 10) +
      scale_color_tq() +
      NoLegend()
      
    ggsave(file.path("plots","umap","ctrl_cytof","clusters_no_legend.pdf"), width = 7.7, height = 5.17, useDingbats = F)

#compartment
    ctrl_cytof %>%
      ggplot(aes(X, Y, color = Region)) +
      geom_point(size=2) +
      theme_void() +
      scale_color_brewer(palette = "Set1") +
      guides(colour = guide_legend(override.aes = list(size=5))) #from url: https://stackoverflow.com/questions/20415963/how-to-increase-the-size-of-points-in-legend-of-ggplot2
    
    ggsave(file.path("plots","umap","ctrl_cytof","regions_with_legend.pdf"), width = 7.7, height = 5.17, useDingbats = F)
    
    ctrl_cytof %>%
      ggplot(aes(X, Y, color = Region)) +
      geom_point(size=2) +
      theme_void() +
      scale_color_brewer(palette = "Set1") +
      NoLegend()
    
    ggsave(file.path("plots","umap","ctrl_cytof","regions_no_legend.pdf"), width = 7.7, height = 5.17, useDingbats = F)

#comparment marimekko
    mosaicGG2(ctrl_cytof, "ClusterID", "Region") +
      scale_fill_brewer(palette = "Set1")
    ggsave(file.path("plots","others","ctrl_cytof","regions_clusters_marimekko.pdf"), width = 7.7, height = 5.17, useDingbats = F)

    hyper_test_n(ctrl_cytof, var1 = "ClusterID", var2 = "Region")    
        
#make heatmap using seurat
    cyt_obj <- t(ctrl_cytof[,c(8:43)]) %>% 
      as.data.frame() 
    
    colnames(cyt_obj) <- gsub("V","cell_",colnames(cyt_obj))
    cyt_obj <- cyt_obj %>% 
      CreateSeuratObject() %>% 
      SCTransform() %>% 
      RunPCA() %>% 
      FindNeighbors() %>% 
      RunUMAP(features=rownames(.)) %>% 
      FindClusters(resolution = seq(.2,2,by=.2))
      
    pdf(file.path("plots","umap","ctrl_cytof","gene_expr_umaps_seurat.pdf"), width = 7.7, height = 5.17, useDingbats = F)
    walk(rownames(cyt_obj), function(x) {
      plt <- FeaturePlot(cyt_obj,x)+
        theme_void() +
        scale_color_viridis_c()
      print(plt)
    })
    dev.off()
    
    #add the metadata from other analysis
    cyt_obj <- AddMetaData(cyt_obj, data.frame(ctrl_cytof[,c(1:7)], row.names = colnames(cyt_obj)))
    cyt_obj2 <- cyt_obj
    Idents(cyt_obj2) <- cyt_obj2$ClusterID
    
    DimPlot(cyt_obj2, label = T)
    
    #find Markers
    all.markers <- cyt_obj2 %>% 
      FindAllMarkers()
    
    top10 <- all.markers %>% 
      group_by(cluster) %>% 
      top_n(10,wt=avg_log2FC)
    cyt_obj2 %>% 
      DoHeatmap(features = unique(top10$gene),
                group.colors = unname(palette_light())) +
      scale_fill_viridis_c(option = "E")
      
    ggsave(file.path("plots","heatmaps","ctrl_cytof","top10_heatmap.pdf"))
    
    #plotheatmap only
    DoHeatmap(cyt_obj2,features = unique(top10$gene), group.bar = F)+
      theme_void() +
      scale_fill_viridis(option = "E")  + 
      NoLegend() 
    
    ggsave(file.path("plots","heatmaps","ctrl_cytof","top10_heatmap_only.png"))
    
    
    #plot protein expression
    walk(c("CD206","CD163","CD169","HLA.DR") , function(x) {
       .df <- data.frame("UMAP_1"=cyt_obj2$X, "UMAP_2"=cyt_obj2$Y,as.data.frame(t(as.matrix(cyt_obj[["SCT"]]@scale.data)))) 
       .df <- .df %>% 
         dplyr::arrange(.df[[`x`]])
        plt <- .df  %>% 
          ggplot(aes(UMAP_1, UMAP_2, color=.data[[`x`]])) +
          geom_point(size=4) +
          scale_color_gradientn("Expr.",colors = colorRampPalette(c("darkblue","lightblue2","yellow","red2"))(100)) +
          theme_void() +
          labs(title = x)
      print(plt)
      ggsave(file.path("plots","umap","ctrl_cytof",paste(x,"expr_umap.pdf", sep = "_")), width = 7.7, height = 5.17, useDingbats = F)
      
      })
    
    #violin plots
    df <- data.frame(ClusterID=cyt_obj2$ClusterID,ctrl_cytof[c(8:43)]) #as.data.frame(t(as.matrix(cyt_obj2[["SCT"]]@scale.data))))
    colnames(df)[2:37] <- colnames(ctrl_cytof)[c(8:43)]
    proteins <- colnames(df)[-1]
    
    cytof_red <- df %>%
      pivot_longer(CD45:CD61 ,"Protein", "Expression")
    
    proteins_clust <- hclust(dist(t(scale(df[,proteins]))))
    colnames(df[,proteins])[proteins_clust$order]
    
    cytof_red$Protein <- factor(cytof_red$Protein, levels = colnames(df[,proteins])[proteins_clust$order])
    
    ggplot(cytof_red, aes(Protein, value+.1, fill=ClusterID)) +
      geom_violin(scale = "width", lwd=0.5) +
      facet_grid(facets = ~ClusterID,
                 drop = TRUE, 
                 scales = "free_y", 
                 space = "free") +
      coord_flip() +
      theme_pubr() +
      labs(y="Expression (A.U.)",x=element_blank()) +
      scale_fill_tq() +
      scale_y_continuous(breaks = c(0,5,10)) +
      theme(
        strip.text.x = element_text(
          size = 12, color = "black", face="bold"
        ),
        legend.position = "None")
    
    ggsave(file.path("plots","others","ctrl_cytof","gene-violin-plots-cytof.pdf"), height = 12, width = 12, useDingbats=F)
    
