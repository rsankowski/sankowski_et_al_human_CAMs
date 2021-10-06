library(RaceID)
library(Seurat)
library(tidyverse)
library(assertthat)
library(FateID)
library(readxl)
library(RColorBrewer)

load(file.path("data", "men_cp_cams_seurat_10x.RData"))
source(file.path("R","functions.R"))

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

date <- Sys.Date()

set.seed(79106)

if (!file.exists(file.path("data","analysis_output_10x_cams", "sc.RData"))) {
  
  sc <- SCseq(cams[["SCT"]]@counts)
  
  #filter data
  sc <- filterdata(sc, 
                   mintotal=100,
                   minnumber = 1,
                   knn=10,
                   minexpr = 1)
  
  # 2.Run Seurat with filtering such that the same cells are retained
  assert_that(length(colnames(cams)) == length(colnames(sc@ndata)))
  
  # 3.Re-initialize RaceID output with Seurat data:
  
  part <- as.numeric(as.character(cams@meta.data$seurat_clusters))
  d <- as.matrix( dist(cams@reductions$pca@cell.embeddings) )
  tsne <- as.data.frame(cams@reductions$umap@cell.embeddings)
  names(part) <- colnames(sc@ndata)
  
  n <- colnames(sc@ndata)
  part <- part[n]
  
  # partition
  sc@cpart <- sc@cluster$kpart <- part
  # distances
  sc@distances <- d[n,n]
  # tsne
  sc@tsne <- tsne[n,]
  rm(d, tsne)
  
  sc@medoids <- compmedoids(sc, sc@cpart)
  
  #reorder the clusters
  idx <- which(order_clusters %in% unique(as.character(sc@cpart)))
  sc@cpart <- factor(sc@cpart, levels = order_clusters)
  
  sc@fcol <- c(colors_pat, colors_many, colors_fig)[-2][idx]
  
  save(sc, file=file.path("data","analysis_output_10x_cams", "sc.RData"))
  
} else {
  load(file.path("data","analysis_output_10x_cams", "sc.RData"))
}


#StemID
if (!file.exists(file.path("data","analysis_output_10x_cams", "ltr.RData"))){ #data/ltr-larger-clusters.RData
  ltr <- Ltree(sc)
  
  #convert clusters in integers
  ltr@sc@cpart <- as.numeric(as.character(ltr@sc@cpart)) +1
  names(ltr@sc@cpart) <- colnames(ltr@sc@ndata)
  
  ltr <- compentropy(ltr)
  ltr <- projcells(ltr,nmode=TRUE,fr=FALSE) #400
  ltr <- projback(ltr,pdishuf=100)
  ltr <- lineagegraph(ltr)
  ltr <- comppvalue(ltr,pthr=0.2)
  
  save(ltr, file = file.path("data","analysis_output_10x_cams", "ltr.RData"))
} else {
  load(file.path("data","analysis_output_10x_cams", "ltr.RData"))
}

#plot entropy umap
.ent <- data.frame(ltr@sc@tsne, Entropy=ltr@entropy)

plot_continuous(param = "Entropy")
ggsave(file.path("plots","umap","10x_cams","entropy.pdf"), width = 11.2, height = 8.11, useDingbats = F)

x <- compscore(ltr,scthr=0.2)
plotgraph(ltr,showCells=FALSE)

pdf(file.path("plots","umap","10x_cams","lineage-graph-tsne-plot.pdf"), width = 11.2, height = 8.11, useDingbats = F)
plotgraph(ltr,showCells=FALSE,cex=4)
dev.off()

x <- compscore(ltr,scthr=0.9)
plotdistanceratio(ltr)
plotspantree(ltr)
plotprojections(ltr)

#lineage tree for moDCs
cams <- c(2,3,1)

#pseudotemporal cams
n <- cellsfromtree(ltr,cams)
x <- getfdata(ltr@sc)

fs  <- filterset(x,n=n$f, minexpr = 2, minnumber = 2)

if (!file.exists("data/s1d-cams.Robj")) {
  s1d <- getsom(fs,nb=1000,alpha=.5)
  save(s1d, file = "data/s1d-cams.Robj")
} else {
  load("data/s1d-cams.Robj")
}
ps  <- procsom(s1d,corthr=.75,minsom=5)
y    <- ltr@sc@cpart[n$f]
fcol <- sc@fcol
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)

pdf(file.path("plots","heatmaps","10x_cams","cams-trajectory-heatmap.pdf"), width = 8.57, height = 5.79, useDingbats = F)
plotheatmap(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

pdf(file.path("plots","heatmaps","10x_cams","outline-cams-trajectory-heatmap.pdf"), width = 8.57, height = 5.79, useDingbats = F)
plotheatmap2(ps$all.z,xpart=y,xcol=fcol,ypart=ps$nodes,xgrid=FALSE,ygrid=TRUE,xlab=FALSE)
dev.off()

png(file.path("plots","heatmaps","10x_cams","map-cams-trajectory-heatmap.png"))
image(t(as.matrix(ps$all.z)), col = rev(colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100)), axes = FALSE, ylim = c(-0.02,1))
dev.off()


#export node genes
modules <- data.frame('Node' = NA, 'Genes' = NA)
for (i in 1: max(ps$nodes)) {
  gene_names <- names(ps$nodes)[ps$nodes == i]
  gene_names <- gsub('_.*', '', gene_names)
  modules2 <- data.frame('Node' = NA, 'Genes' = gene_names)
  modules2$Node <- rep(as.character(i), nrow(modules2))
  modules <- rbind(na.omit(modules), modules2)
}

write_csv(modules, file.path("data","analysis_output_10x_cams", "nodes-stemid-vector-cams.csv"))
