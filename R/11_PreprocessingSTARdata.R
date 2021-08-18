library(readr)
library(Matrix)
##### NOTE: transfer all the Solo output folders in the Star_SoloOut folder!  #####
listF <- list.dirs(file.path("data","Star_SoloOut"),recursive = F)

Allcounts <- c()
for(i in listF){
  GeneDir <- file.path(  i, "Gene", "raw")
  genes <- readr::read_delim(file.path(GeneDir, "features.tsv"), col_names = FALSE, delim = "\t")
  cells <- readr::read_delim(file.path(GeneDir, "barcodes.tsv"), col_names = FALSE, delim = "\t")
  n_genes <- nrow(genes)
  n_cells <- nrow(cells)
  matr <- readr::read_delim(file.path(GeneDir, "matrix.mtx"), skip = 3, delim = " ", col_names = FALSE)
  counts <- Matrix::sparseMatrix(
    i = matr$X1,
    j = matr$X2,
    x = matr$X3,
    dims = c(n_genes, n_cells)
  )
  rownames(counts)  <-genes$X2
  colnames(counts) <- paste0(i, "_",cells$X1)
  Allcounts <- cbind(Allcounts[genes$X2, ],counts )
}

Allspliced <- c()
Allunspliced <- c()
Allambiguous <- c()
for(i in listF){
  VeloDir <- file.path(i, "Velocyto", "raw")
  genes <- readr::read_delim(file.path(VeloDir, "features.tsv"), col_names = FALSE, delim = "\t")
  cells <- readr::read_delim(file.path(VeloDir, "barcodes.tsv"), col_names = FALSE, delim = "\t")
  n_genes <- nrow(genes)
  n_cells <- nrow(cells)
  spl <- readr::read_delim(file.path(VeloDir, "spliced.mtx"), skip = 3, delim = " ", col_names = FALSE)
  spliced <- Matrix::sparseMatrix(
    i = spl$X1,
    j = spl$X2,
    x = spl$X3,
    dims = c(n_genes, n_cells)
  )
  unspl <- readr::read_delim(file.path(VeloDir, "unspliced.mtx"), skip = 3, delim = " ", col_names = FALSE)
  unspliced <- Matrix::sparseMatrix(
    i = unspl$X1,
    j = unspl$X2,
    x = unspl$X3,
    dims = c(n_genes, n_cells)
  )
  amb <- readr::read_delim(file.path(VeloDir, "ambiguous.mtx"), skip = 3, delim = " ", col_names = FALSE)
  ambiguous <- Matrix::sparseMatrix(
    i = amb$X1,
    j = amb$X2,
    x = amb$X3,
    dims = c(n_genes, n_cells)
  )
  rownames(spliced) <- rownames(unspliced) <- rownames(ambiguous) <-genes$X2
  colnames(spliced) <- colnames(unspliced) <- colnames(ambiguous)<- paste0(i, "_",cells$X1)
  Allspliced <- cbind(Allspliced[genes$X2, ],spliced )
  Allunspliced <- cbind(Allunspliced[genes$X2, ],unspliced )
  Allambiguous <- cbind(Allambiguous[genes$X2, ],ambiguous )
}



BrainData_List <- list(counts = Allcounts,
                           spliced = Allspliced,
                           unspliced = Allunspliced,
                           ambiguous = Allambiguous)

#adjust cell IDs and filter for cells in cams
load(file.path("data", "men_cp_seurat_10x.RData"))

BrainData_List <- map(BrainData_List, function(x) {
  colnames(x) <- case_when(
    grepl("data/Star_SoloOut/3-men_Solo.out_", colnames(x)) ~ paste0(colnames(x),"-1_2"),
    T ~ paste0(colnames(x),"-1_1")
  )
  colnames(x) <- gsub("(data/Star_SoloOut/3-men_Solo.out_|data/Star_SoloOut/4-cp_Solo.out_)","",colnames(x))
  x1 <- x[,colnames(all)]
  x1
})

#check that cell IDs exist in the 10x Seurat object
assert_that(sum(colnames(BrainData_List$unspliced) %in% colnames(all)) == dim(all)[2])

save(BrainData_List, file=file.path("data","10x_velopcyto_output","BrainData_List.RData"))
