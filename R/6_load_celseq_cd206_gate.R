# adjusted from Josip Herman
library(readr)
library(data.table)
library(tidyverse)
library(tools)
library(assertthat)
library(SingleCellExperiment)
library(scDblFinder)
library(Matrix)

# Load necessary files
gene2iso <- read.csv(file.path("data","wgEncodeGencodeBasicVM9_clean_genes2groups.tsv"),sep="\t",header=FALSE)
ercc     <- read.csv(file.path("data","ERCC_Controls_Analysis_length_ext.txt"),sep="\t",header=TRUE)
geneid_ercc    <- data.frame("GENEID" = c(as.character(gene2iso$V1), as.character(ercc$Re.sort.ID)), stringsAsFactors = F )

# Specify directories
file_paths <- file.path("data","counts_CD206_gate")

# Search for prdata files in the specified directories and extract files by a given extension
files_list            <- list.files(file_paths, pattern = ".coutt.csv")

# Get complete paths to prdata files
prdata_file_paths        <- file.path(file_paths, files_list)
names(prdata_file_paths) <- lapply(strsplit(prdata_file_paths,split="/"), function(x) { sub(".coutt.csv","",x[length(x)]) } )

# Calculate md5sums to check for duplicated
md5sums      <- md5sum(prdata_file_paths)

# Check for duplicated data
assert_that(sum(duplicated(md5sums)) == 0)

# Check for duplicated names
assert_that(sum(duplicated(prdata_file_paths)) == 0)

####
#### LOADING
####
# Loading data using lapply
data_list   <- map(prdata_file_paths, function(x) {fread(x, header= T)} )

#add control microglia
file_paths            <- file.path("data","counts_nn")
files_list            <- list.files(file_paths, pattern = ".coutt.csv")

# Get complete paths to prdata files
prdata_file_paths        <- file.path(file_paths, files_list)
names(prdata_file_paths) <- lapply(strsplit(prdata_file_paths,split="/"), function(x) { sub(".coutt.csv","",x[length(x)]) } )

# Calculate md5sums to check for duplicated
md5sums      <- md5sum(prdata_file_paths)

# Check for duplicated data
assert_that(sum(duplicated(md5sums)) == 0)

# Check for duplicated names
assert_that(sum(duplicated(prdata_file_paths)) == 0)

if (F) {
####
#### LOADING
####
# Loading data using lapply
data_list2   <- map(prdata_file_paths, function(x) {fread(x, header= T)} )

#anonymize the names of list
names(data_list2) <- paste("nn", names(data_list2), sep="_")

#add control microglia
data_list <- c(data_list, data_list2)
}

# Add dataset name prefix to prdata columns, Merge with remaining gene names
for (d in names(data_list)) { 
  colnames(data_list[[d]]) <- c("GENEID", paste(d, "_",1:192,sep="" ))
}



# Cbind list of data.tables and removing the GENEID column from data.tables
data_list_cbind <- purrr::reduce(data_list, full_join, by = "GENEID")
data_list_cbind <- data_list_cbind %>% full_join(geneid_ercc)
data_list_cbind[is.na(data_list_cbind)]   <- 0

# Measure time
Sys.time() - start.time

#get the data frame and adjust row names
prdata <- as.data.frame(data_list_cbind)
rownames(prdata) <- prdata$GENEID
prdata$GENEID    <- NULL

#Remove cells with less than 200 umis and genes that were nor detected
prdata <- prdata[names(rowSums(prdata)[rowSums(prdata)>0]) ,names(colSums(prdata)[colSums(prdata)>=200])]

#adjust rownames
genes <- which(!duplicated(gsub("_.*", "", rownames(prdata))))
prdata <- prdata[genes,]
rownames(prdata) <- gsub("_.*", "", rownames(prdata))

#remove doublets
sce <- SingleCellExperiment(assays=list(counts=Matrix(as.matrix(prdata), sparse = T)))
sce <- scDblFinder(sce)
prdata <- prdata[,sce@colData$scDblFinder.class == "singlet"]

save(prdata, file = file.path("data","prdata_celseq_adult_cd206_gate_only.RData"))
saveRDS(prdata, file.path("data","prdata_celseq_adult_cd206_gate_only.rds"))
