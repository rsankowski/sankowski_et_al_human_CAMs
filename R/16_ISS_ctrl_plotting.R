#giotto microlgia only
library(Giotto)
library(tidyverse)
library(viridis)
library(tidyquant)
library(data.table)

source(file.path("R","functions.R"))

#set wd
my_working_dir <- file.path(getwd(),"data")
python_path = "/Users/romansankowski/opt/anaconda3/bin/python"

cell_type_order <- c('CAMs', 'microglia', 'cDCs_monocytes', 'astrocytes', 'oligodendrocytes','OPC', 'ependymal_cells', 'Neurons','VLMC', 'vascular_endothelial_cells', "VSMC", "vascular_endothelial_and_mural_cells",'pericytes')

b1hi_ctrl_cells <- read.csv(file.path("data","iss_ctrl_gbm","b1hi_ctrl_cells_final_prior08.csv"), row.names = 1) %>% 
  mutate(cell_type = factor(cell_type, levels = cell_type_order))
b1hi_ctrl_reads <- read_csv(file.path("data","iss_ctrl_gbm","segmentation_counts_filt_prior08.csv")) %>%
  as.data.frame()

#introduce rownames
rownames(b1hi_ctrl_reads) <- b1hi_ctrl_reads[,1]
reads <- b1hi_ctrl_reads[,-1]

#heatmap
a <- CreateSeuratObject(reads, min.features = 5) %>% 
  SCTransform() %>% 
  RunPCA() %>% 
  RunUMAP(features=rownames(.)) %>% 
  FindNeighbors() %>% 
  FindClusters()

a$cell_type <- b1hi_ctrl_cells[colnames(a), "cell_type"]

Idents(a) <- a$cell_type

if (!file.exists(file.path("data", "markers_all.RData"))) {
  markers_all <- FindAllMarkers(a,only.pos = T)          
  save(markers_all, file = file.path("data", "markers_all.RData"))
  write_csv(markers_all, file.path("data", "markers_all.csv"))
} else {
  load(file.path("data", "markers_all.RData"))
}

top10 <- markers_all %>% 
  #remove non annotated genes
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

heat <- DoHeatmap(a,features = top10$gene, group.colors = unname(palette_light())) 
heat + 
  scale_fill_viridis(option = "F")

ggsave(file.path("plots", "heatmaps", "iss_ctrl","clusters_all_heatmap.pdf"), useDingbats=F, height=16, width=20)

DoHeatmap(a,features = top10$gene, group.bar = FALSE)  + 
  scale_fill_viridis(option = "F") +
  theme_void() +
  NoLegend()

ggsave(file.path("plots", "heatmaps", "iss_ctrl", "clusters_all_heatmap_hetamap_only.png"), height=16, width=20)

#other code

#exclude cells with less than 5 reads 
reads <- reads[, which(colSums(b1hi_ctrl_reads[,-1])>5)]

df_all <- reads %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var="cell") %>%
  left_join(b1hi_ctrl_cells[, c("cell", "cell_type")]) %>%
  pivot_longer(MOG:PLSCR1, names_to = "Gene", values_to="Expression")
colnames(df_all)[1] <- c("ID")

ref <- read_csv2(file.path("data","iss_ctrl_gbm","Cell_type_markers_cartana_CAMs.csv")) %>%
  mutate(Cell_type = gsub("(vascular_and_leptomeningeal_cells_|vascular_smooth_muscle_cells_)", "", Cell_type),
         Cell_type = as.factor(Cell_type))

ref$Cell_type <- factor(ref$Cell_type, levels = c(cell_type_order, "multiple_cell_types", "Interferon_genes"))

reads2 <- split(reads[ref$Gene,], ref$Cell_type)
reads2 <- lapply(reads2, function(x) {
  a <- rowSums(x)
  a <- rownames_to_column(as.data.frame(a), var="gene")
  a <- reorder(a$gene, desc(a$a))
  a
})
reads2 <- unlist(reads2)

ref %>% group_by()
df_all$cell_type <- factor(df_all$cell_type, levels = c(cell_type_order, "multiple types of glia", "Interferon_genes"))
df_all$Gene <- factor(df_all$Gene, levels = rev(levels(reads2)))


df_all %>%
  ggplot(aes(x=ID, y=Gene, fill=log(Expression+.1))) +
  geom_tile() +
  facet_grid(facets = ~cell_type, 
             drop = TRUE, 
             #space = "free", 
             scales = "free", 
             #switch = "x",
             space = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.spacing.x=unit(0.05, "lines"),
        panel.spacing.y=unit(1, "lines")) +
  labs(x=element_blank()) +
  scale_fill_viridis(option="B")

ggsave("plots/single_cell_heatmap_block1_ctrl_hi_thresh.pdf")

df_all %>%
  ggplot(aes(x=ID, y=Gene, fill=log(Expression+.1))) +
  geom_tile() +
  facet_grid(facets = ~cell_type, 
             drop = TRUE, 
             #space = "free", 
             scales = "free", 
             #switch = "x",
             space = "free_x") +
  theme_void() +
  theme(axis.text.x = element_blank(),
        panel.spacing.x=unit(0.05, "lines"),
        panel.spacing.y=unit(1, "lines")) +
  scale_fill_viridis(guide=F, option = "B")

ggsave("plots/heatmap_clustering_CAMs_block1_ctrl_hi_thresh_image_only.png")

#find variable genes
ids <- df_all[! duplicated(df_all[, c("ID", "cell_type")]),c("ID", "cell_type")]
vars <- split(as.data.frame(t(reads)), ids$cell_type)
vars2 <- lapply(vars[1:13], function(x) {
  x1 <- as.data.frame(x)
  a <- apply(x, 2, sd)
  b <- names(a)[a>.5]
  b <- b[b %in% ref$Gene[ref$defining == "yes"]]
  x1[,b]
  
})

var <- apply(reads[,unique(df_all$ID[df_all$cell_type == "microglia"])], 1, sd)
names(var)[var>.5]

#cluster data 
#CAMs
#determine optimal number of cluster
set.seed(31)
# function to compute total within-cluster sum of squares
df_scaled <- scale(vars2[[1]])
fviz_nbclust(df_scaled, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")
fviz_nbclust(df_scaled, kmeans, method = "silhouette", k.max = 24) + theme_minimal() + ggtitle("The Silhouette Plot")

a <- kmeans(df_scaled, 6)$cluster
clust <- data.frame(ID=names(a), cluster=a)

df_all[df_all$Gene %in% colnames(vars2[[1]]),] %>%
  left_join(clust) %>%
  na.omit() %>%
  group_by(cluster, Gene) %>%
  summarise(Expression= mean(Expression)) %>%
  ggplot(aes(x=factor(cluster), y= Gene, fill=log(Expression+.1))) +
  #stat_summary(fun=mean, geom = "tile") +
  geom_tile(height=1, width=1) +
  scale_fill_viridis_c(option = "A") +
  theme_minimal() +
  theme(panel.grid = element_blank()) + 
  coord_fixed(ratio = 1) +
  labs(x="Cluster", y="Gene")

ggsave("plots/heatmap_clustering_CAMs_block1_ctrl_hi_thresh.pdf")

mt <- df_all[df_all$Gene %in% colnames(vars2[[1]]),] %>%
  left_join(clust) %>%
  na.omit() %>%
  group_by(cluster, Gene) %>%
  summarise(Expression= mean(Expression)) %>%
  pivot_wider(names_from = "Gene", values_from="Expression") %>%
  as.matrix()

rownames(mt) <- mt[,1]

pheatmap::pheatmap(t(mt[,-1]), scale="row")

#micr
micr_scaled <- scale(vars2[[2]])
fviz_nbclust(micr_scaled, kmeans, method = "wss", k.max = 24) + theme_minimal() + ggtitle("the Elbow Method")
fviz_nbclust(micr_scaled, kmeans, method = "silhouette", k.max = 24) + theme_minimal() + ggtitle("The Silhouette Plot")

a <- kmeans(micr_scaled, 5)$cluster

clust <- data.frame(ID=names(a), cluster=a)

df_all[df_all$Gene %in% colnames(vars2[[2]]),] %>%
  left_join(clust) %>%
  na.omit() %>%
  group_by(cluster, Gene) %>%
  summarise(Expression= mean(Expression)) %>%
  ggplot(aes(x=factor(cluster), y= Gene, fill=log(Expression+.1))) +
  #stat_summary(fun=mean, geom = "tile") +
  geom_tile(height=1, width=1) +
  scale_fill_viridis_c(option = "A") +
  theme_minimal() +
  theme(panel.grid = element_blank()) + 
  coord_fixed(ratio = 1) +
  labs(x="Cluster", y="Gene")

ggsave("plots/heatmap_clustering_microglia_block1_ctrl_hi_thresh.pdf")

#get representative images
#filter the dataset for cells with more than 5 reads
cells2 <- b1hi_ctrl_cells %>%   
  filter(n_transcripts>=5)

#split neighboring cells into groups of 10

if (!file.exists(file.path("data","iss_ctrl_gbm","ctrl_clusters.RData"))) {
  set.seed(79106)
  .cells <- cells2[, c("x","y")]
  clusters <- kmeans(.cells, centers = nrow(.cells)/12)$cluster
  
  save(clusters, file = file.path("data","iss_ctrl_gbm","ctrl_clusters.RData"))
} else {
  load(file.path("data","iss_ctrl_gbm","ctrl_clusters.RData"))
}

cells2 <- cells2 %>% 
  mutate(cluster = clusters)

  cells3 <- cells2 %>% 
    distinct(cluster, cell_type, .keep_all = T) %>%
    group_by(cluster) %>% 
    mutate(cell_type2=paste(cell_type,collapse = ",")) %>% 
    filter(grepl("(CAMs.*microglia.*astrocytes.*VLMC|microglia.*CAMs)", cell_type2))
  
  cells4 <- cells3 %>% 
    distinct(cluster, .keep_all = T) %>%
    rowwise %>%
    #adjusted from https://stackoverflow.com/questions/54277205/search-for-multiple-values-in-a-column-in-r
    mutate(SUM = str_split(cell_type2, ",", simplify = T) %>% 
             map( ~ str_count(cell_type_order[c(1,2,3,4,9,11)], .)) %>%
             unlist %>% sum) %>% 
    filter(SUM >3)

  cells2_sub <- cells2 %>% 
    filter(cluster %in% cells4$cluster)

  #plot sample clusters
  walk(cells4$cluster[1:10], function(i) {
      plt <- cells2 %>% 
        dplyr::filter(cluster %in% as.character(i)) %>% 
        ggplot(aes(x,y,color=cell_type)) +
        geom_point()
      print(plt)
    })
  
#load dapi image  
  dapi <- fread("/Users/romansankowski/Documents/other_analyses/20201125_cartana_spatial/data/Image_data_pairs/block_1/DAPI_Block 1_greyscale_CARTANA8730401_20201102.txt")
  
  #load raw reads
  raw_reads <- read.csv("/Users/romansankowski/Documents/other_analyses/20201125_cartana_spatial/data/Image_data_pairs/block_1/reads_Block 1_HighThreshold_CARTANA8730401_20201102.csv")

  #define markers and colors
  set.seed(79106)
  marker_col <- expand.grid(c(0:6,8,15:18), unname(palette_dark())) %>% 
    sample_n(length(unique(raw_reads$gene))) %>% 
    transmute(pch=Var1,
              col=as.character(Var2),
              gene = unique(raw_reads$gene)) 
  
  raw_reads <- raw_reads %>% 
    left_join(marker_col)
  
  walk(1:nrow(cells4), 
       function(i) {
         .cluster <- cells4[i, "cluster"]
         .cells2 <- cells2[cells2$cluster == as.numeric(.cluster),]
         .raw_reads <- raw_reads[which(raw_reads$X > (min(.cells2$x) - 20) & raw_reads$X < (max(.cells2$x) + 20) & raw_reads$Y > (min(.cells2$y) - 20) & raw_reads$Y < (max(.cells2$y) + 20)), ]
         .pixel <- dapi[as.numeric(rownames(dapi)) > (min(.cells2$x) - 50) & as.numeric(rownames(dapi)) < (max(.cells2$x) + 50)] %>% 
               as.data.frame() %>% 
               select(colnames(dapi)[as.numeric(gsub("V", "",colnames(dapi))) > (min(.cells2$y) - 50) & as.numeric(gsub("V", "",colnames(dapi))) < (max(.cells2$y) + 50)]) 
          .pixel$x = rownames(dapi)[as.numeric(rownames(dapi)) > (min(.cells2$x) - 50) & as.numeric(rownames(dapi)) < (max(.cells2$x) + 50)]
          .pixel <- .pixel %>% 
            pivot_longer(-x, names_to = "y", values_to = "intensity") %>% 
            mutate(y=as.numeric(y), y=as.numeric(gsub("V","",.$y)),
                   x=as.numeric(x))
          plt <- .raw_reads %>% ggplot(aes(X,Y,color=factor(gene), shape=factor(gene))) +
                geom_raster(data = .pixel, aes(x,y,fill=intensity), inherit.aes = F) + 
                #scale_fill_gradientn(colours = c("black","skyblue")) +
                #geom_point(data = .cells2, aes(x, y), inherit.aes = F, color="navy", size=5) +
                geom_point(size=3) +
                #adjusted from url https://stackoverflow.com/questions/12410908/combine-legends-for-color-and-shape-into-a-single-legend
                scale_colour_manual(name = "Gene",
                                    labels = marker_col$gene,
                                    values = as.character(marker_col$col)) +   
                scale_shape_manual(name = "Gene",
                                   labels = marker_col$gene,
                                   values = marker_col$pch) + 
                theme_void() +
                geom_text(data=.cells2, aes(x,y,label=cell_type), inherit.aes = F) + #, color = "white"
                labs(title = i)
          print(plt)
          
       })
  
  .raw_reads %>% ggplot(aes(X,Y,color=factor(gene), shape=factor(gene))) +
    geom_raster(data = .pixel, aes(x,y,fill=intensity), inherit.aes = F) + 
    geom_point(size=3) +
    #adjusted from url https://stackoverflow.com/questions/12410908/combine-legends-for-color-and-shape-into-a-single-legend
    scale_colour_manual(name = "Gene",
                        labels = marker_col$gene,
                        values = as.character(marker_col$col)) +   
    scale_shape_manual(name = "Gene",
                       labels = marker_col$gene,
                       values = marker_col$pch) + 
    theme_void() +
    geom_text(data=.cells2, aes(x,y,label=cell_type), inherit.aes = F) 
