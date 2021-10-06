#(devtools)
# install_github("velocyto-team/velocyto.R")
# http://velocyto.org/
library("velocyto.R")
library(Seurat)
require(scales)
library(pagoda2)

#load functions
source(file.path("R","functions.R"))

load(file.path("data","10x_velocyto_output","BrainData_List.RData"))
load(file.path("data", "men_cp_cams_seurat_10x.RData"))

# get embedding
embUMAP <- cams[["umap"]]@cell.embeddings
embPCA <- cams[["pca"]]@cell.embeddings

### define cluster order
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

### set colors
clusters <- cams$seurat_clusters
# unique(clusters)
palette <- c(colors_pat, colors_many, colors_fig)[-2][1:length(levels(cams))]
names(palette) <- levels(clusters)
cell.colors <- palette[clusters]
names(cell.colors) <- colnames(cams)


#### subset loom file
# exonic read (spliced) expression matrix
NoDupl <- !duplicated(rownames(BrainData_List$counts))
emat <- BrainData_List$spliced[, colnames(cams)]
# intronic read (unspliced) expression matrix
nmat <- BrainData_List$unspliced[, colnames(cams)]
# spanning read (intron+exon) expression matrix
# smat <- ldat$spanning;

hist(log10(rowSums(emat)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')
hist(log10(colSums(emat)),col='wheat',xlab='cell size')


# filter expression matrices based on some minimum max-cluster averages
emat2 <- filter.genes.by.cluster.expression(emat,cell.colors,min.max.cluster.average = 0.5)
nmat2 <- filter.genes.by.cluster.expression(nmat,cell.colors,min.max.cluster.average = 0.05)
# look at the resulting gene set
length(intersect(rownames(emat2),rownames(nmat2)))


#run velocyto
cell.dist <- as.dist(1-armaCor(t(cams@reductions$pca@cell.embeddings)))

fit.quantile <- 0.02
rvel.cd <- gene.relative.velocity.estimates(emat2,nmat2,deltaT=1,kCells=20, 
                                            min.nmat.smat.correlation = 0.05, min.nmat.emat.correlation = 0.05,
                                            min.nmat.emat.slope = 0.05,
                                            cell.dist=cell.dist,fit.quantile=fit.quantile)
# to save time, follow this:
# https://github.com/velocyto-team/velocyto.R/issues/15
x <- show.velocity.on.embedding.cor(embUMAP,rvel.cd,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),
                               cex=0.8,arrow.scale=1.2,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=23,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0.1)



show.velocity.on.embedding.cor(embUMAP,rvel.cd,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.5),
                               cex=0.8,arrow.scale=5,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=40,
                               arrow.lwd=1,do.par=F,cell.border.alpha = 0.1, cc = x$cc)


show.velocity.on.embedding.cor(embUMAP,rvel.cd,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.8),
                               cex=0.8,arrow.scale=2,arrow.lwd=1,fixed.arrow.length=T,cc = x$cc) +
  theme_void()

pdf(file.path("plots","umap","10x_cams","velocity_umap.pdf"), useDingbats = F, width = 11.2, height = 8.11)
show.velocity.on.embedding.cor(embUMAP,rvel.cd,n=300,scale='sqrt',cell.colors=ac(cell.colors,alpha=0.2),
                               cex=4,arrow.scale=2,arrow.lwd=1.5,fixed.arrow.length=T,cc = x$cc)
dev.off()

# Plot spliced and unspliced UMAPs for specific genes
gene <- "CD163"
GenePlot <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 20,kGenes=1,
                                 fit.quantile=fit.quantile,cell.emb=embUMAP,
                                 cell.colors=cell.colors,cell.dist=cell.dist,
                                 show.gene=gene,old.fit=rvel.cd,do.par=T)



