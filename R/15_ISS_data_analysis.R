library(Giotto)
library(tidyverse)
library(smfishHmrf)
library(tidyquant)

source(file.path("R","functions.R"))

#set wd
my_working_dir <- file.path(getwd(),"data", "iss_ctrl_gbm")
python_path = "/Users/romansankowski/opt/anaconda3/bin/python"

#assign cell types
ref <- read_csv2(file.path("data","iss_ctrl_gbm","Cell_type_markers_cartana_CAMs.csv")) %>%
  filter(defining == "yes" ) %>%
  mutate(Cell_type = gsub("(vascular_and_leptomeningeal_cells_|vascular_smooth_muscle_cells_)", "", Cell_type),
         Cell_type = as.factor(Cell_type))

#extract expression matrix
b1hi_ctrl_reads <- read.delim(file.path("data","iss_ctrl_gbm","segmentation_counts.tsv"), sep = "\t") %>%
  as.data.frame()

rownames(b1hi_ctrl_reads) <- b1hi_ctrl_reads$gene
b1hi_ctrl_reads <- b1hi_ctrl_reads[,-1]
colnames(b1hi_ctrl_reads) <- gsub("X","cell", colnames(b1hi_ctrl_reads))

#split data frame according to cell type
df4 <- b1hi_ctrl_reads[ref$Gene,]

df4 <- split(df4, ref$Cell_type, drop=T)
df5 <- lapply(df4, colSums)
df5 <- bind_rows(df5, .id="Cell_type")
df5 <- as.data.frame(df5)         
rownames(df5) <- df5$Cell_type
df5 <- as.data.frame(t(df5[,-1]))

cell_types <- as.factor(apply(df5, 1, which.max))
levels(cell_types) <- levels(ref$Cell_type)

save(cell_types, file = file.path("data","iss_ctrl_gbm","cell_types_hi_thresh_prior_08.RData"))

#load cell stats
b1hi_ctrl_cells <- read_csv(file.path("data","iss_ctrl_gbm","segmentation_cell_stats.csv")) %>%
  na.omit()
b1hi_ctrl_cells$cell <- paste0("cell",as.character(b1hi_ctrl_cells$cell))
b1hi_ctrl_cells$cell_type <- cell_types[b1hi_ctrl_cells$cell]

#inspect data
b1hi_ctrl_cells %>%
  ggplot(aes(x,y,color=factor(cluster))) +
  geom_point(size=.5)

#cell type
b1hi_ctrl_cells %>%
  ggplot(aes(x,y,color=factor(cell_type))) +
  geom_point(size=.5)

#plot blood vessel cells
walk(unique(b1hi_ctrl_cells$cell_type), function(x) {
  plt <- b1hi_ctrl_cells %>%
    filter(cell_type %in% x) %>%
    ggplot(aes(x,y,color=factor(cell_type))) +
    geom_point(size=.5)
  
  print(plt)
})

#assign sampleID to cells
set.seed(79106)
if (!file.exists(file.path("data","iss_ctrl_gbm","sampleID_b1hi_ctrl_prior08.Rdata"))) {
  sampleID <- kmeans(b1hi_ctrl_cells[, 2:3], centers = 60, iter.max = 20)$cluster
  save(sampleID, file = file.path("data","iss_ctrl_gbm","sampleID_b1hi_ctrl_prior08.Rdata"))
}else{
  load(file.path("data","iss_ctrl_gbm","sampleID_b1hi_ctrl_prior08.Rdata"))
}

b1hi_ctrl_cells$sampleID <- sampleID

#set cluster medians
b1hi_ctrl_cells_med <- b1hi_ctrl_cells %>%
  group_by(sampleID) %>% 
  summarize(x=median(x),
            y=median(y))
#plot sample IDs
b1hi_ctrl_cells %>%
  ggplot(aes(x,y,color=factor(sampleID))) +
  geom_point(size=.5) +
  geom_text(data = b1hi_ctrl_cells_med, aes(x,y,label=sampleID), color="black")


b1hi_ctrl_cells$sampleID <- case_when(
  b1hi_ctrl_cells$sampleID %in% c("16", "10", "55", "32", "12", "38", "18", "4", "36", "15", "18", "60", "1", "37", "17") ~ "ctrl1",
  b1hi_ctrl_cells$sampleID %in% c("30", "34", "2", "24", "39", "41","59","23","43","48","51","19","26","25","3", "52","49") ~ "GBM1",
  b1hi_ctrl_cells$sampleID %in% c("8", "47", "31", "57", "53", "13", "5", "22") ~ "GBM2",
  TRUE ~ "ctrl2"
)

b1hi_ctrl_cells$diagnosis <- gsub("[1-9]", "", b1hi_ctrl_cells$sampleID)
#inspect data
b1hi_ctrl_cells %>%
  ggplot(aes(x,y,color=factor(diagnosis))) +
  geom_point(size=.5)

b1hi_ctrl_cells %>%
  ggplot(aes(x,y,color=factor(sampleID))) +
  geom_point(size=.5)

#extract controls
b1hi_ctrl_cells <- b1hi_ctrl_cells %>%
  filter(diagnosis == "ctrl")

#filter the counts data
b1hi_ctrl_reads[, b1hi_ctrl_cells$cell ] %>%
  write.csv(file.path("data","iss_ctrl_gbm","segmentation_counts_filt_prior08.csv"))

write_csv(b1hi_ctrl_cells,file.path("data","iss_ctrl_gbm","b1hi_ctrl_cells_final_prior08.csv"))

if (!file.exists(file.path("data","iss_ctrl_gbm","b1hi_ctrl_coord_prior08.csv"))) {
  # <- reformat the coordinate file
  b1hi_ctrl_coord <- b1hi_ctrl_cells[,1:3]
  colnames(b1hi_ctrl_coord) <- c("ID", "X", "Y")
  write_csv(b1hi_ctrl_coord, file.path("data","iss_ctrl_gbm","b1hi_ctrl_coord_prior08.csv"))
} else {
  b1hi_ctrl_coord <- read_csv(file.path("data","iss_ctrl_gbm","b1hi_ctrl_coord_prior08.csv"))
}

#run giotto
#tutorial: http://spatialgiotto.rc.fas.harvard.edu/giotto.seqfish.html
# 1. (optional) set Giotto instructions
instrs = createGiottoInstructions(save_plot = TRUE,show_plot = FALSE,save_dir = my_working_dir,python_path = python_path)
# 2. create giotto object from provided paths ####
expr_path = fs::path(my_working_dir, file.path("segmentation_counts_filt_prior08.csv"))
loc_path = fs::path(my_working_dir, file.path("b1hi_ctrl_coord_prior08.csv"))

#create giotto object
iss1hi <- createGiottoObject(raw_exprs = expr_path, spatial_locs = loc_path, instructions = instrs)

#add metadata
iss1hi <- addCellMetadata(iss1hi,new_metadata = b1hi_ctrl_cells,by_column = T,column_cell_ID = 'cell')
cell_metadata = pDataDT(iss1hi)

iss1hi <- filterGiotto(gobject = iss1hi,expression_threshold = 1,gene_det_in_min_cells = 5,min_det_genes_per_cell = 5, expression_values = c('raw'),verbose = T)
## normalize
iss1hi <- normalizeGiotto(gobject = iss1hi, scalefactor = 6000, verbose = T)
## add gene & cell statistics
iss1hi <- addStatistics(gobject = iss1hi)
## adjust expression matrix for technical or known variables
iss1hi <- adjustGiottoMatrix(gobject = iss1hi, expression_values = c('normalized'),batch_columns = NULL, covariate_columns = c('nr_genes', 'total_expr'),return_gobject = TRUE,update_slot = c('custom'))

## visualize
spatPlot(gobject = iss1hi, save_param = list(save_name = '2_spatplot'))

#highly variable genes
iss1hi <- calculateHVG(gobject = iss1hi, method = 'cov_loess', difference_in_cov = 0.001, save_param = list(save_name = '3_a_HVGplot', base_height = 5, base_width = 5))
## select genes based on HVG and gene statistics, both found in gene metadata
gene_metadata = fDataDT(iss1hi)
featgenes = gene_metadata[hvg == 'yes' & perc_cells > 4 & mean_expr_det > 0.5]$gene_ID

#run PCA UMAP etc
# runPCA: normal and recommended usage, set center = T:
iss1hi <- runPCA(gobject = iss1hi, genes_to_use = rownames(b1hi_ctrl_reads), scale_unit = T, center = T, ncp=100, maxit=3000)
screePlot(iss1hi, save_param = list(save_name = '3_b_screeplot'))
plotPCA(gobject = iss1hi,save_param = list(save_name = '3_c_PCA_reduction'))

## run UMAP and tSNE on PCA space (default)
iss1hi <- runUMAP(iss1hi, dimensions_to_use = 1:10, n_threads = 10, genes_to_use = rownames(b1hi_ctrl_reads)) #c("MS4A7", "CD163", "LYVE1", "CSF1R", "SOX9", "TYROBP", "C1QA", "P2RY12", "CX3CR1", "AQP4", "TMEM119","HLA-DRA", "COL1A1", "FCER1A", "TMEM163", "COL23A1", "IFI44", "RBFOX3","SIGLEC1","CTSS", "CCR2", "S100A9","TRIM22", "SLC2A5","SOX2", "DDX60L", "TGFBI","ITGAM", "MRC1",  "CLDN5")
plotUMAP(gobject = iss1hi,save_param = list(save_name = '3_d_UMAP_reduction'))

ucoords <- data.frame(dim1=iss1hi@dimension_reduction$cells$umap$umap$coordinates[,1],
                      dim2=iss1hi@dimension_reduction$cells$umap$umap$coordinates[,2])
b <- ucoords %>%
  bind_cols(data.frame(n_transcripts = iss1hi@cell_metadata$n_transcripts)) %>%
  ggplot(aes(dim1, dim2, color=log(n_transcripts))) +
  geom_point(size=0.5) +
  scale_color_viridis_c()

b

#nearest neighbor
iss1hi <- createNearestNetwork(gobject = iss1hi, dimensions_to_use = 1:10, k = 15)

#leiden clustering 
iss1hi <- doLeidenCluster(gobject = iss1hi, resolution = .2, n_iterations = 1000)

#plot
plotUMAP(gobject = iss1hi,cell_color = 'leiden_clus', show_NN_network = F, point_size = 0.5,save_param = list(save_name = '4_a_UMAP_leiden'))
plotUMAP(gobject = iss1hi,cell_color = 'cell_type', show_NN_network = F, point_size = 0.5,save_param = list(save_name = '4_a_UMAP_leiden')) +
  facet_wrap(~ cell_type, ncol=4)

#find marker genes
gini_markers_subclusters = findMarkers_one_vs_all(gobject = iss1hi,method = 'gini', expression_values = 'normalized',cluster_column = 'leiden_clus', min_genes = 20, min_expr_gini_score = 0.5, min_det_gini_score = 0.5)
topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']
# violinplot
violinPlot(iss1hi, genes = unique(topgenes_gini$genes), cluster_column = 'leiden_clus',strip_text = 8, strip_position = 'right', cluster_custom_order = unique(topgenes_gini$cluster),save_param = c(save_name = '6_a_violinplot_gini', base_width = 5, base_height = 10))
# cluster heatmap
topgenes_gini2 = gini_markers_subclusters[, head(.SD, 6), by = 'cluster']
plotMetaDataHeatmap(iss1hi, selected_genes = unique(topgenes_gini2$genes), custom_gene_order = unique(topgenes_gini2$genes),custom_cluster_order = unique(topgenes_gini2$cluster),metadata_cols = c('leiden_clus'), x_text_size = 10, y_text_size = 10,save_param = c(save_name = '6_b_metaheatmap_gini'))
ggsave(file.path("plots", "heatmaps", "iss_ctrl", "metadata_heatmap.pdf"), useDingbats=F)

#plot cell type
cell_types2 <- cell_types[iss1hi@cell_metadata$cell_ID]
names(cell_types2) <- iss1hi@cell_metadata$cell_ID
cell_type_order <- c('CAMs', 'microglia', 'astrocytes', 'oligodendrocytes','OPC', 'ependymal_cells', 'Neurons','VLMC', 'vascular_endothelial_cells', "VSMC", "vascular_endothelial_and_mural_cells",'pericytes', 'cDCs_monocytes')

iss1hi <- annotateGiotto(gobject = iss1hi, annotation_vector = cell_types2, cluster_column = 'cell_ID', name = 'cell_types')

assertthat::are_equal(iss1hi@cell_metadata$cell_type, iss1hi@cell_metadata$cell_types)

#correct
gini_markers_subclusters[, cell_types := cell_types[cluster] ]
gini_markers_subclusters[, cell_types := factor(cell_types, cell_type_order)]
data.table::setorder(gini_markers_subclusters, cell_types)

#create spatial grid
iss1hi <- createSpatialGrid(gobject = iss1hi,sdimx_stepsize = 500,sdimy_stepsize = 500,minimum_padding = 50)
spatPlot(gobject = iss1hi, show_grid = T, point_size = 1.5,save_param = c(save_name = '8_a_grid'))

#delaunay network
plotStatDelaunayNetwork(gobject = iss1hi, maximum_distance = 400, save_plot = F)
iss1hi = createSpatialNetwork(gobject = iss1hi, minimum_k = 2, maximum_distance_delaunay = 400)

#9.2. Create spatial networks based on k and/or distance from centroid
iss1hi <- createSpatialNetwork(gobject = iss1hi, method = 'kNN', k = 5, name = 'spatial_network')
iss1hi <- createSpatialNetwork(gobject = iss1hi, method = 'kNN', k = 10, name = 'large_network')
iss1hi <- createSpatialNetwork(gobject = iss1hi, method = 'kNN', k = 100,
                               maximum_distance_knn = 200, minimum_k = 2, name = 'distance_network')

#visualize network data
cell_metadata = pDataDT(iss1hi)
#sample ctrl2
field1_ids = cell_metadata[sampleID == "ctrl2"]$cell_ID
subiss1hi = subsetGiotto(iss1hi, cell_ids = field1_ids)
spatPlot(gobject = subiss1hi, show_network = T,network_color = 'blue', spatial_network_name = 'Delaunay_network',point_size = 2.5, cell_color = 'cell_types', save_param = c(save_name = '9_a_spatial_network_delaunay', base_height = 6))
spatPlot(gobject = subiss1hi, show_network = T,network_color = 'blue', spatial_network_name = 'spatial_network',point_size = 2.5, cell_color = 'cell_types',save_param = c(save_name = '9_b_spatial_network_k3', base_height = 6))
spatPlot(gobject = subiss1hi, show_network = T,network_color = 'blue', spatial_network_name = 'large_network',point_size = 2.5, cell_color = 'cell_types',save_param = c(save_name = '9_c_spatial_network_k10', base_height = 6))
spatPlot(gobject = subiss1hi, show_network = T,network_color = 'blue', spatial_network_name = 'distance_network',point_size = 2.5, cell_color = 'cell_types',save_param = c(save_name = '9_d_spatial_network_dist', base_height = 6))

#spatial genes
km_spatialgenes = binSpect(iss1hi)
spatGenePlot(iss1hi, expression_values = 'scaled', genes = km_spatialgenes[1:4]$genes,point_shape = 'border', point_border_stroke = 0.1,show_network = F, network_color = 'lightgrey', point_size = 2.5, cow_n_col = 2,save_param = list(save_name = '10_a_spatialgenes_km'))

#10.2. Spatial genes co-expression modules
ext_spatial_genes = km_spatialgenes[adj.p.value<.05]$genes

spat_cor_netw_DT = detectSpatialCorGenes(iss1hi, method = 'network', spatial_network_name = 'Delaunay_network',subset_genes = ext_spatial_genes)
spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT, name = 'spat_netw_clus', k = 8)
heatmSpatialCorGenes(iss1hi, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',save_param = c(save_name = '10_b_spatialcoexpression_heatmap',base_height = 6, base_width = 8, units = 'cm'), heatmap_legend_param = list(title = NULL), show_row_names = T)

pdf(file.path("plots", "heatmaps", "iss_ctrl", "block1_ctrl_hithresh_heatmSpatialCorGenes.pdf"), height = 7, width = 9, useDingbats = F)
heatmSpatialCorGenes(iss1hi, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',save_param = c(save_name = '10_b_spatialcoexpression_heatmap',base_height = 8, base_width = 8, units = 'cm'), heatmap_legend_param = list(title = NULL), show_row_names = T)
dev.off()

netw_ranks = rankSpatialCorGroups(iss1hi, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',save_param = c(save_name = '10_c_spatialcoexpression_rank',base_height = 3, base_width = 5))
top_netw_spat_cluster = showSpatialCorGenes(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',selected_clusters = 6, show_top_genes = 1)
cluster_genes_DT = showSpatialCorGenes(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus', show_top_genes = 1)
cluster_genes = cluster_genes_DT$clus; names(cluster_genes) = cluster_genes_DT$gene_ID
iss1hi = createMetagenes(iss1hi, gene_clusters = cluster_genes, name = 'cluster_metagene')
spatCellPlot(iss1hi,spat_enr_names = 'cluster_metagene',cell_annotation_values = netw_ranks$clusters,point_size = 1.5, cow_n_col = 3,save_param = c(save_name = '10_d_spatialcoexpression_metagenes',base_width = 11, base_height = 6))

#11. HMRF spatial domains

hmrf_folder = fs::path(my_working_dir,'11_HMRF/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
my_spatial_genes = km_spatialgenes[adj.p.value<.05]$genes
# do HMRF with different betas
HMRF_spatial_genes = doHMRF(gobject = iss1hi, expression_values = 'scaled',spatial_genes = my_spatial_genes,spatial_network_name = 'Delaunay_network',k = 9,betas = c(28,2,3), output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_top100_k9_scaled'))
## view results of HMRF
for(i in seq(28, 32, by = 2)) {
  viewHMRFresults2D(gobject = iss1hi,HMRFoutput = HMRF_spatial_genes,k = 9, betas_to_view = i,point_size = 2)
}
## add HMRF of interest to giotto object
iss1hi = addHMRF(gobject = iss1hi,HMRFoutput = HMRF_spatial_genes,k = 9, betas_to_add = c(28),hmrf_name = 'HMRF_2')
## visualize
spatPlot(gobject = iss1hi, cell_color = 'HMRF_2_k9_b.28', point_size = 3, coord_fix_ratio = 1, save_param = c(save_name = '11_HMRF_2_k9_b.28', base_height = 3, base_width = 9, save_format = 'png'))

#12. Cell neighborhood
#12.1. cell-type/cell-type interactions

cell_proximities = cellProximityEnrichment(gobject = iss1hi,cluster_column = 'cell_types',spatial_network_name = 'Delaunay_network',adjust_method = 'BH')
## barplot
cellProximityBarplot(gobject = iss1hi,CPscore = cell_proximities, min_orig_ints = 4, min_sim_ints = 4, save_param = c(save_name = '12_a_barplot_cell_cell_enrichment'))
## heatmap
cellProximityHeatmap(gobject = iss1hi, CPscore = cell_proximities, order_cell_types = T, scale = T,color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'),save_param = c(save_name = '12_b_heatmap_cell_cell_enrichment', unit = 'in'))
pdf(file.path("plots", "heatmaps", "iss_ctrl", "block1_ctrl_hithresh_cellProximityHeatmap.pdf"), useDingbats = F)
cellProximityHeatmap(gobject = iss1hi, CPscore = cell_proximities, order_cell_types = T, scale = T,color_breaks = c(-1.5, 0, 1.5), color_names = c('blue', 'white', 'red'),save_param = c(save_name = '12_b_heatmap_cell_cell_enrichment', unit = 'in'))
dev.off()
## network
cellProximityNetwork(gobject = iss1hi, CPscore = cell_proximities, remove_self_edges = T,only_show_enrichment_edges = T,save_param = c(save_name = '12_c_network_cell_cell_enrichment'))
pdf(file.path("plots", "others", "iss_ctrl", "block1_ctrl_hithresh_cellProximityNetwork.pdf"), useDingbats = F)
cellProximityNetwork(gobject = iss1hi, CPscore = cell_proximities, remove_self_edges = T,only_show_enrichment_edges = T,save_param = c(save_name = '12_c_network_cell_cell_enrichment'))
dev.off()

## network with self-edges
cellProximityNetwork(gobject = iss1hi, CPscore = cell_proximities,remove_self_edges = F, self_loop_strength = 0.3,only_show_enrichment_edges = F,rescale_edge_weights = T,node_size = 8,edge_weight_range_depletion = c(1, 2),edge_weight_range_enrichment = c(2,5),save_param = c(save_name = '12_d_network_cell_cell_enrichment_self',base_height = 5, base_width = 5, save_format = 'pdf'))
## visualization of specific cell types
# Option 1
spec_interaction = "CAMs--microglia"
cellProximitySpatPlot2D(gobject = iss1hi,interaction_name = spec_interaction,show_network = T,cluster_column = 'cell_types',cell_color = 'cell_types',cell_color_code = c(astrocytes = 'lightblue', Olig = 'red'),point_size_select = 4, point_size_other = 2,save_param = c(save_name = '12_e_cell_cell_enrichment_selected'))
# Option 2: create additional metadata
iss1hi = addCellIntMetadata(iss1hi, spatial_network = 'spatial_network',cluster_column = 'cell_types',cell_interaction = spec_interaction,name = 'astro_olig_ints')
spatPlot(iss1hi, cell_color = 'astro_olig_ints',select_cell_groups =  c('other_astrocytes', 'other_Olig', 'select_astrocytes', 'select_Olig'),legend_symbol_size = 3, save_param = c(save_name = '12_f_cell_cell_enrichment_sel_vs_not'))

#12.2. Cell neighborhood: interaction changed genes

## select top 25th highest expressing genes
gene_metadata = fDataDT(iss1hi)
plot(gene_metadata$nr_cells, gene_metadata$mean_expr)
plot(gene_metadata$nr_cells, gene_metadata$mean_expr_det)
quantile(gene_metadata$mean_expr_det)
high_expressed_genes = gene_metadata[mean_expr_det > 1.31]$gene_ID
## identify genes that are associated with proximity to other cell types
CPGscoresHighGenes =  findICG(gobject = iss1hi,selected_genes = high_expressed_genes,spatial_network_name = 'Delaunay_network',cluster_column = 'cell_types',diff_test = 'permutation',adjust_method = 'fdr',nr_permutations = 2000, do_parallel = T, cores = 2)
## visualize all genes
plotCellProximityGenes(iss1hi, cpgObject = CPGscoresHighGenes, method = 'dotplot', save_param = c(save_name = '13_a_CPG_dotplot', base_width = 5, base_height = 5))
## filter genes
CPGscoresFilt = filterICG(CPGscoresHighGenes, min_log2_fc = 1, min_fdr = 0.1)

## visualize subset of interaction changed genes (ICGs)
ICG_genes = c("P2RY12","SPP1","CSF1R","MFGE8","S100A9","COL23A1","MYH11","COL1A1") #unique(CPGscoresFilt$CPGscores$genes) "CAVIN1","COL15A1" 
ICG_genes_types = as.factor(c("microglia", "astrocytes", "cDCs_monocytes","VSMC","VLMC"))
names(ICG_genes) = ICG_genes_types
plotICG(gobject = iss1hi,cpgObject = CPGscoresHighGenes,source_type = 'CAMs',source_markers = c("LYVE1", "HLA-DRA","C3"),ICG_genes = ICG_genes,save_param = c(save_name = '13_b_ICG_barplot'),
        cell_color_code = unname(palette_light()[c(4,3,2,9,11)]))

ggsave(file.path("plots", "others", "iss_ctrl", "cell_interaction_barplot_CAMs_micr_astrocytes.pdf"), useDingbats=F, width = 12, height = 4)

#export data
save(iss1hi, file =file.path("data", "iss_ctrl_gbm", "iss1_all_block_1_hi_thresh.RData"))

