#index data celseq2
library(data.table)
library(tidyverse)
library(viridis)
library(ggpubr)
library(Seurat)
library(FSA)

source(file.path("R","functions.R"))

#load data
load(file.path("data","SCORE-integration.RData"))
metadata <- bind_cols(all[[]], UMAP_1= all@reductions$NetUMAP@cell.embeddings[,1], UMAP_2= all@reductions$NetUMAP@cell.embeddings[,2])
metadata$cell_ID <- rownames(metadata)

#set order for clusters
order_clusters <- data.frame(seurat_clusters= all@meta.data[,"seurat_clusters"], row.names = rownames(all@meta.data)) %>%
        bind_cols(as.data.frame(t(all[["SCT"]]@scale.data))) %>%
        group_by(seurat_clusters) %>%
        summarize_all(.funs=mean) %>%
        as.data.frame()

rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder clusters
metadata$seurat_clusters <- factor(metadata$seurat_clusters, levels = rev(order_clusters))

#Load datasets(
if (!file.exists(file.path("data","Index_anon","index_all_ctrl.RData"))) {
        #plate layout
        bcd_plate <- read_csv(file.path("data","384-Well-plate-layout.csv"))
        
        index_list <- list()
        
        #Pat 24
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat24_Men'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_1_', bcd_plate2$Barcode), paste0(sample_name, '_2_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat24_Men_1.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat24", Region_ind = "Other")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat 12
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat12_PVM'
        bcd_plate2$cell_ID <- paste0(sample_name, '_6_', bcd_plate2$Barcode)
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat12_PVM_1_gm.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat12", Region_ind = "GM")
        index2 <- read_csv("data/Index_anon/Pat12_PVM_1_wm.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat12", Region_ind = "WM")
        index <- bind_rows(index1, index2)
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat4_dura
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat4_Dura_8'
        bcd_plate2$cell_ID <- paste0(sample_name, "_",bcd_plate2$Barcode)
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat4_Dura_1.csv")  %>%
                add_column(Compartment_ind = "Dura", Diagnosis_ind = "Ctrl", Pat_ID = "Pat4", Region_ind = "Other")
        index2 <- bind_rows(read_csv("data/Index_anon/Pat4_Dura_2.csv")) %>%
                add_column(Compartment_ind = "Dura", Diagnosis_ind = "Ctrl", Pat_ID = "Pat4", Region_ind = "Other")
        index3 <- bind_rows(read_csv("data/Index_anon/Pat4_Dura_3.csv")) %>%
                add_column(Compartment_ind = "Dura", Diagnosis_ind = "Ctrl", Pat_ID = "Pat4", Region_ind = "Other")
        index4 <- bind_rows(read_csv("data/Index_anon/Pat4_Dura_4.csv")) %>%
                add_column(Compartment_ind = "Dura", Diagnosis_ind = "Ctrl", Pat_ID = "Pat4", Region_ind = "Other")
        
        index <- index1 %>%
                bind_rows(index2) %>%
                bind_rows(index3) %>%
                bind_rows(index4)
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat2_Men
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat2_Men'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_1_3_', bcd_plate2$Barcode), paste0(sample_name, '_2_4_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat2_Men.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat2", Region_ind = "Other")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat8_Dura_Men
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat8_Dura_Men'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_3_', bcd_plate2$Barcode), paste0(sample_name, '_4_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- tibble(read_csv("data/Index_anon/Pat8_Dura_Men_1.csv"),Compartment_ind = "Dura", Diagnosis_ind = "Ctrl", Pat_ID = "Pat25", Region_ind = "Other") %>%
                bind_rows(tibble(read_csv("data/Index_anon/Pat8_Dura_Men_2.csv"),Compartment_ind = "Dura", Diagnosis_ind = "Ctrl", Pat_ID = "Pat25", Region_ind = "Other")) %>%
                bind_rows(tibble(read_csv("data/Index_anon/Pat8_Dura_Men_3.csv"),Compartment_ind = "Dura", Diagnosis_ind = "Ctrl", Pat_ID = "Pat26", Region_ind = "Other")) %>%
                bind_rows(tibble(read_csv("data/Index_anon/Pat8_Dura_Men_4.csv"),Compartment_ind = "MenM",Diagnosis_ind = "Ctrl", Pat_ID = "Pat9", Region_ind = "Other")) 
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat14_CP
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat14_CP_8'
        bcd_plate2$cell_ID <- paste0(sample_name, '_', bcd_plate2$Barcode)
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat14_CP.csv")  %>%
                add_column(Compartment_ind = "CP", Diagnosis_ind = "Ctrl", Pat_ID = "Pat14", Region_ind = "Other")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat14_Men
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat14_Men'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_9_', bcd_plate2$Barcode), paste0(sample_name, '_2_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat14_Men.csv")  %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat14", Region_ind = "Other")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat14_PVM
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat14_PVM'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_5_', bcd_plate2$Barcode), paste0(sample_name, '_2_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat14_PVM_GM.csv")  %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat14", Region_ind = "GM")
        index <- read_csv("data/Index_anon/Pat14_PVM_WM.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat14", Region_ind = "WM") %>%
                bind_rows(index1)
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat21_GM_WM_PVM
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat21_GM_WM_PVM'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_5_', bcd_plate2$Barcode), paste0(sample_name, '_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat21_GM_WM_PVM_1.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat21", Region_ind = "GM")
        
        index2 <- read_csv("data/Index_anon/Pat21_GM_WM_PVM_2.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat21", Region_ind = "GM")
        
        index <- read_csv("data/Index_anon/Pat21_GM_WM_PVM_3.csv") %>%
                filter(`Sort Index X`<13) %>% 
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat21", Region_ind = "GM") %>% 
                bind_rows(index2) %>% 
                bind_rows(index1)
        
        
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat21_Men
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat21_Men_2'
        bcd_plate2$cell_ID <- paste0(sample_name, '_', bcd_plate2$Barcode)
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat21_Men.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat21", Region_ind = "Other")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat15_Dura
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat15_Dura_7'
        bcd_plate2$cell_ID <- paste0(sample_name, '_', bcd_plate2$Barcode)
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat15_Dura.csv") %>%
                add_column(Compartment_ind = "Dura", Diagnosis_ind = "Ctrl", Pat_ID = "Pat15", Region_ind = "Other")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat3_Men_9
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat3_Men_9'
        bcd_plate2$cell_ID <- paste0(sample_name, '_', bcd_plate2$Barcode)
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat3_Men_9.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat3", Region_ind = "Other")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat5_Dura
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat5_Dura_1'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_1_', bcd_plate2$Barcode), paste0(sample_name, '_2_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat5_Dura_1_1.csv") %>%
                bind_rows(read_csv("data/Index_anon/Pat5_Dura_1_2.csv")) %>%
                bind_rows(read_csv("data/Index_anon/Pat5_Dura_1_3.csv")) %>%
                add_column(Compartment_ind = "Dura", Diagnosis_ind = "Ctrl", Pat_ID = "Pat5", Region_ind = "Other")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #PVM1
        bcd_plate2 <- bcd_plate
        sample_name <- 'PVM_7'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_2_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/PVM1_1.csv")%>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat45", Region_ind = "GM")
        index2 <- read_csv("data/Index_anon/PVM1_2.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat45", Region_ind = "WM")
        
        index3 <- read_csv("data/Index_anon/PVM1_3.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat31", Region_ind = "GM")
        
        index4 <- read_csv("data/Index_anon/PVM1_4.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat31", Region_ind = "WM")
        
        index5 <- read_csv("data/Index_anon/PVM1_5.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat3", Region_ind = "GM")
        
        index <- read_csv("data/Index_anon/PVM1_6.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat3", Region_ind = "WM") %>%
                bind_rows(index1) %>%
                bind_rows(index2) %>%
                bind_rows(index3) %>%
                bind_rows(index4) %>%
                bind_rows(index5)
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat18_PVM_1
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat18_PVM_4'
        bcd_plate2$cell_ID <- paste0(sample_name, '_', bcd_plate2$Barcode)
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat18_PVM_1.csv")  %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat18", Region_ind = "GM")
        
        index <- read_csv("data/Index_anon/Pat18_PVM_2.csv")  %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat18", Region_ind = "WM") %>%
                bind_rows(index1)
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2[bcd_plate2$Y>12,] %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat20_CP
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat20_CP_4'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_1_', bcd_plate2$Barcode), paste0(sample_name, '_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat20_CP.csv")  %>%
                add_column(Compartment_ind = "CP", Diagnosis_ind = "Ctrl", Pat_ID = "Pat20", Region_ind = "Other")
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat20_GM_PVM
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat20_GM_PVM_6'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_1_', bcd_plate2$Barcode), paste0(sample_name, '_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat20_GM_PVM_1.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat20", Region_ind = "GM")
        index <- read_csv("data/Index_anon/Pat21_GM_WM_PVM_3.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat21", Region_ind = "WM") %>%
                bind_rows(index1)
        
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat20_GM_WM_PVM
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat20_GM_WM_PVM_9'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_2_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat20_GM_WM_PVM_1.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat20", Region_ind = "WM")
        index <- read_csv("data/Index_anon/Pat20_GM_WM_PVM_2.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat20", Region_ind = "GM")  %>%
                bind_rows(index1)
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat20_Men
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat20_Men_3'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_2_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat20_Men.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat20", Region_ind = "Other")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat1_Men
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat1_Men'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_1_3_', bcd_plate2$Barcode), paste0(sample_name, '_2_4_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat1_Men.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat1", Region_ind = "Other")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #CP5
        bcd_plate2 <- bcd_plate
        sample_name <- 'CP_5'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_2_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/CP5.csv") %>%
                add_column(Compartment_ind = "CP", Diagnosis_ind = "Ctrl", Pat_ID = "Pat36", Region_ind = "Other")
        index <- read_csv("data/Index_anon/CP5_2.csv")  %>%
                add_column(Compartment_ind = "CP", Diagnosis_ind = "Ctrl", Pat_ID = "Pat44", Region_ind = "Other") %>% 
                bind_rows(index1)
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat16_GM_PVM
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat16_GM_PVM_3'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_2_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat16_GM_PVM_1.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat16", Region_ind = "GM")
        index <- read_csv("data/Index_anon/Pat16_GM_PVM_2.csv")  %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat16", Region_ind = "WM") %>%
                bind_rows(index1)
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat16_Pat18_Men_1.csv
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat16_Pat18_Men_4'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_1_', bcd_plate2$Barcode), paste0(sample_name, '_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat16_Pat18_Men_1.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat16", Region_ind = "Other")
        index <- read_csv("data/Index_anon/Pat16_Pat18_Men_2.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat18", Region_ind = "Other") %>%
                bind_rows(index1)
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat16_17_CP
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat16_Pat17_CP_3'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_2_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat16_Pat17_CP_1.csv") %>%
                add_column(Compartment_ind = "CP", Diagnosis_ind = "Ctrl", Pat_ID = "Pat16", Region_ind = "Other")
        
        index2 <- read_csv("data/Index_anon/Pat16_Pat17_CP_2.csv") %>%
                add_column(Compartment_ind = "CP", Diagnosis_ind = "Ctrl", Pat_ID = "Pat36", Region_ind = "Other")
        
        index <- read_csv("data/Index_anon/Pat16_Pat17_CP_3.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat18", Region_ind = "GM") %>%
                bind_rows(index1) %>%
                bind_rows(index2)
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat9_GM_PVM
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat9_GM_PVM_5'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_2_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat9_GM_PVM.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat9", Region_ind = "GM")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #Pat9_WM_PVM
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat9_WM_PVM_6'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_1_', bcd_plate2$Barcode), paste0(sample_name, '_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat9_WM_PVM.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat9", Region_ind = "WM")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #GBM Samples "Pat19_GBM_PVM"
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat19_GBM_PVM'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_1_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat19_GBM_PVM.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "GBM", Pat_ID = "Pat19", Region_ind = "WM")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #"Pat11_GBM_Micr_1"
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat11_GBM_Micr'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_1_', bcd_plate2$Barcode), paste0(sample_name, '_2_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat11_GBM_Micr.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "GBM", Pat_ID = "Pat19", Region_ind = "WM")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        
        
        
        #Pat6_10_GBM_Men
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat6_10_GBM_Men'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_1_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat6_10_GBM_Men.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "GBM", Pat_ID = "Pat6", Region_ind = "Other")
        
        index <- read_csv("data/Index_anon/Pat6_10_GBM_Men_2.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "GBM", Pat_ID = "Pat11", Region_ind = "Other") %>%
                bind_rows(index1)
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        
        
        
        #Pat13_GBM_PVM 
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat13_GBM_PVM'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_1_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat13_GBM_PVM.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "GBM", Pat_ID = "Pat13", Region_ind = "Other")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        
        
        
        #"Pat22_GBM_Micr 
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat22_GBM_Micr'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_1_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat22_GBM_Micr.csv") %>%
                add_column(Compartment_ind = "Micr", Diagnosis_ind = "GBM", Pat_ID = "Pat22", Region_ind = "WM")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        
        #"Pat22_GBM_PVM
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat22_GBM_PVM'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_1_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index <- read_csv("data/Index_anon/Pat22_GBM_PVM.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "GBM", Pat_ID = "Pat22", Region_ind = "WM")
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        
        #'Pat7_Pat23_GBM_PVM', 
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat7_Pat23_GBM_PVM'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_1_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat7_Pat23_GBM_PVM_1.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "GBM", Pat_ID = "Pat23", Region_ind = "WM")
        
        index <- read_csv("data/Index_anon/Pat7_Pat23_GBM_PVM_2.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "GBM", Pat_ID = "Pat7", Region_ind = "WM") %>%
                bind_rows(index1)
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        
        #'#'Pat10_11_GBM_PVM
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat10_11_GBM_PVM'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_1_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat10_11_GBM_PVM_1.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "GBM", Pat_ID = "Pat10", Region_ind = "WM")
        
        index <- read_csv("data/Index_anon/Pat10_11_GBM_PVM_2.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "GBM", Pat_ID = "Pat11", Region_ind = "WM") %>%
                bind_rows(index1)
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #'#Pat8_GBM_Men' 
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat8_GBM_Men'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_1_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat8_GBM_Men_1.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat27", Region_ind = "Other")
        
        index2 <- read_csv("data/Index_anon/Pat8_GBM_Men_2.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat11", Region_ind = "Other")
        
        index3 <- read_csv("data/Index_anon/Pat8_GBM_Men_3.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat29", Region_ind = "Other")
        
        index4 <- read_csv("data/Index_anon/Pat8_GBM_Men_4.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "GBM", Pat_ID = "Pat6", Region_ind = "Other")
        
        index5 <- read_csv("data/Index_anon/Pat8_GBM_Men_5.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "GBM", Pat_ID = "Pat6", Region_ind = "Other")
        
        index <- read_csv("data/Index_anon/Pat8_GBM_Men_6.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "GBM", Pat_ID = "Pat11", Region_ind = "Other") %>%
                bind_rows(index1) %>%
                bind_rows(index2) %>%
                bind_rows(index3) %>%
                bind_rows(index4) %>%
                bind_rows(index5)
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #'#Pat8_GBM_PVM_Men'
        bcd_plate2 <- bcd_plate
        sample_name <- 'Pat8_GBM_PVM_Men'
        bcd_plate2$cell_ID <- ifelse(bcd_plate2$Y <13, paste0(sample_name, '_', bcd_plate2$Barcode), paste0(sample_name, '_1_', bcd_plate2$Barcode))
        bcd_plate2 <- bcd_plate2 %>% left_join(metadata[,!grepl("(Pat_ID|_ind)$", colnames(metadata))])
        
        #load index data
        index1 <- read_csv("data/Index_anon/Pat8_GBM_PVM_Men_1.csv") %>%
                add_column(Compartment_ind = "MenM", Diagnosis_ind = "GBM", Pat_ID = "Pat11", Region_ind = "Other")
        
        index2 <- read_csv("data/Index_anon/Pat8_GBM_PVM_Men_2.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat27", Region_ind = "GM")
        
        index3 <- read_csv("data/Index_anon/Pat8_GBM_PVM_Men_3.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat29", Region_ind = "GM")
        
        index4 <- read_csv("data/Index_anon/Pat8_GBM_PVM_Men_4.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat11", Region_ind = "GM")
        
        index5 <- read_csv("data/Index_anon/Pat8_GBM_PVM_Men_5.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "Ctrl", Pat_ID = "Pat28", Region_ind = "GM")
        
        index <- read_csv("data/Index_anon/Pat8_GBM_PVM_Men_6.csv") %>%
                add_column(Compartment_ind = "PVM", Diagnosis_ind = "GBM", Pat_ID = "Pat10", Region_ind = "Other") %>%
                bind_rows(index1) %>%
                bind_rows(index2) %>%
                bind_rows(index3) %>%
                bind_rows(index4) %>%
                bind_rows(index5)
        
        colnames(index)[1:2] <- c('Y', 'X')
        
        #Merge datasets
        index_full <- bcd_plate2 %>% left_join(index)
        colnames(index_full)
        
        #Adjust the colnames
        a <- colnames(index_full)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub(".*? (.+)", "\\1", a)
        a <- gsub("-", "_", a)
        
        a
        colnames(index_full) <- a
        
        index_list[[sample_name]] <- index_full[index_full$cell_ID %in% rownames(metadata),]
        
        #bind datasets
        index_all <- index_list %>%
                bind_rows(.id="Sample") %>%
                na.omit
        
        index_all_facs <- index_all[, grepl('(Area|Compartment_ind|Diagnosis_ind|Pat_ID|Region_ind|UMAP_1|UMAP_2|seurat_clusters)', colnames(index_all))]
        index_all <- bind_cols(index_all[,1:16],as.data.frame(index_all_facs))
        colnames(index_all) <- gsub("_Area", "", colnames(index_all))
        colnames(index_all)[23:24] <- c("FSC", "SSC")
        
        #adjust MHCII outlier cells
        index_all$MHCII[index_all$MHCII > .3] <- max(index_all$MHCII[index_all$MHCII < .3])
        
        #assign rownames
        index_all <- as.data.frame(index_all)
        rownames(index_all) <- index_all$cell_ID
        
        save(index_all, file = file.path("data","Index_anon","index_all_ctrl.RData"))
} else {
        load(file.path("data","Index_anon","index_all_ctrl.RData"))
}

#plot index
my_cols <- toupper(rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')))

for (i in colnames(index_all)[20:29]) {
        pl <- plot_index(gene=c(i), point_size = 4) + scale_color_gradientn(colors=my_cols) + theme(text=element_text(size=20))
        print(pl)
        ggsave(file.path("plots","umap","celseq_cd206",paste0(i,"-logscale.pdf")), useDingbats=F)
}

for (i in colnames(index_all)[20:29]) {
        pl <- plot_index(gene=c(i), point_size = 4, log=F)  + scale_color_gradientn(colors=my_cols) + theme(text=element_text(size=20))
        print(pl)
        ggsave(file.path("plots","umap","celseq_cd206",paste0(i,"-linear-scale.pdf")), useDingbats=F)
}

#violin plots of index data
index_all <- index_all[,-c(23:24,27)]

index_long <- index_all %>%
        pivot_longer(CD11b:CD206 ,"Protein", "Expression")

proteins <- hclust(dist(t(scale(index_all[,c(21:26)]))))
colnames(index_all[,21:26])[proteins$order]

index_long$Protein <- factor(index_long$Protein, levels = colnames(index_all[,21:26])[proteins$order])

ggplot(index_long, aes(factor(seurat_clusters, levels=rev(levels(metadata$seurat_clusters))), value, fill=seurat_clusters)) +
        geom_violin(scale = "width", lwd=0.5, draw_quantiles = 0.5) +
        facet_grid(facets = ~Protein,
                   drop = TRUE, 
                   #space = "free", 
                   scales = "free_x", 
                   #switch = "x",
                   #space = "free_x"
        ) +
        coord_flip() +
        theme_pubclean() +
        labs(y="Expression (A.U.)",x=element_blank()) +
        scale_fill_manual(values = c(colors_many), guide="none") + 
        #scale_y_continuous(limits=c(0,9), breaks = c(0,5)) +
        theme(
                strip.text.x = element_text(
                        size = 12, color = "black", face="bold"
                ))

ggsave(file.path("plots","others","celseq_cd206","gene-violin-plots-linear-scale.pdf"), height = 7, width = 7,useDingbats=F)

#violin plots log scale
ggplot(index_long, aes(factor(seurat_clusters, levels=rev(levels(metadata$seurat_clusters))), log(value+0.1), fill=seurat_clusters)) +
        geom_violin(scale = "width", lwd=0.5, draw_quantiles = 0.5) +
        facet_grid(facets = ~Protein,
                   drop = TRUE, 
                   #space = "free", 
                   scales = "free_x", 
                   #switch = "x",
                   #space = "free_x"
        ) +
        coord_flip() +
        #theme_pubr() +
        theme_pubclean() +
        labs(y="Expression (A.U.)",x=element_blank()) +
        scale_fill_manual(values = c(colors_many), guide="none") + 
        #scale_y_continuous(limits=c(0,9), breaks = c(0,5)) +
        theme(
                strip.text.x = element_text(
                        size = 12, color = "black", face="bold"
                )) 

ggsave(file.path("plots","others","celseq_cd206","gene-violin-plots-log-scale.pdf"), height = 7, width = 7,useDingbats=F)


#statistical analysis of clusters
#define comparisons
my_comparisons <- list(c("1", "0"), c("1", "3"), c("1", "4"), c("1", "5"), c("1", c("2", "7", "8", "9", "10","6")))

stats <- list()
for (i in colnames(index_all)[21:26]) {
        
        df <- index_all[,c("seurat_clusters",`i`)]
        colnames(df)[2] <- "Protein"
        pl <-  ggplot(df, aes(factor(seurat_clusters, levels=levels(metadata$seurat_clusters)), log(Protein+0.1), fill=seurat_clusters)) +
                geom_violin(scale = "width", lwd=0.5, draw_quantiles = 0.5) +
                #coord_flip() +
                #theme_pubr() +
                theme_pubclean() +
                labs(y="log(Expression (A.U.))",x=element_blank(), title = i) +
                scale_fill_manual(values =  colors_many, guide="none") + 
                #scale_y_continuous(limits=c(0,9), breaks = c(0,5)) +
                theme(
                        strip.text.x = element_text(
                                size = 12, color = "black", face="bold"
                        )) +
                stat_compare_means(comparisons = my_comparisons, size=7)
        
        print(pl)
        ggsave(file.path("plots","others","celseq_cd206",paste0(i,"-gene-violin-plots-log-scale.pdf")), useDingbats=F)
        
        mod <- dunnTest(log(Protein) ~ seurat_clusters, data = df, method = "bh")$res #tidy(kruskal.test(gene ~ Diagnosis, data = panela_down3[panela_down3$Cluster == i,]))
        #filter(P.adj<0.05)
        stats[[i]] <- mod
}

stats <- bind_rows(stats, .id="Protein") %>%
        mutate(P.adj = p.adjust(P.unadj, method = "hochberg")) %>%
        filter(P.adj<.05) %>%
        write_csv(file.path("data","Index_anon","stat-testing-protein-clusters.csv"))

for (i in colnames(index_all)[21:26]) {
        
      
        df <- index_all[,c("seurat_clusters",`i`)]
        colnames(df)[2] <- "Protein"
        pl <-  ggplot(df, aes(factor(seurat_clusters, levels=levels(metadata$seurat_clusters)), Protein, fill=seurat_clusters)) +
                geom_violin(scale = "width", lwd=0.5, draw_quantiles = 0.5) +
                #coord_flip() +
                #theme_pubr() +
                theme_pubclean() +
                labs(y="Expression (A.U.)",x=element_blank(), title = i) +
                scale_fill_manual(values = colors_many, guide=F) + 
                #scale_y_continuous(limits=c(0,9), breaks = c(0,5)) +
                theme(
                        strip.text.x = element_text(
                                size = 12, color = "black", face="bold"
                        ))
        print(pl)
        ggsave(file.path("plots","others","celseq_cd206",paste0(i,"-gene-violin-plots-linear-scale.pdf")), useDingbats=F)
        
}

for (i in colnames(index_all)[21:26]) {
        
        df <- index_all[,c("seurat_clusters","celltype",`i`)]
        colnames(df)[3] <- "Protein"
        pl <-  ggplot(df, aes(factor(celltype), log(Protein+0.1), fill=celltype)) +
                geom_violin(scale = "width", lwd=0.5, draw_quantiles = 0.5) +
                #coord_flip() +
                #theme_pubr() +
                theme_pubclean() +
                labs(y="log(Expression (A.U.))",x=element_blank(), title = i) +
                scale_fill_manual(values =  colors_many, guide="none") + 
                #scale_y_continuous(limits=c(0,9), breaks = c(0,5)) +
                theme(
                        strip.text.x = element_text(
                                size = 12, color = "black", face="bold"
                        )) +
                stat_compare_means(comparisons = my_comparisons, size=7)
        
        print(pl)
        ggsave(file.path("plots","others","celseq_cd206",paste0(i,"-gene-violin-plots-log-scale.pdf")), useDingbats=F)
        
        mod <- dunnTest(log(Protein) ~ seurat_clusters, data = df, method = "bh")$res #tidy(kruskal.test(gene ~ Diagnosis, data = panela_down3[panela_down3$Cluster == i,]))
        #filter(P.adj<0.05)
        stats[[i]] <- mod
}

stats <- bind_rows(stats, .id="Protein") %>%
        mutate(P.adj = p.adjust(P.unadj, method = "hochberg")) %>%
        filter(P.adj<.05) %>%
        write_csv(file.path("data","Index_anon","stat-testing-protein-clusters.csv"))

for (i in colnames(index_all)[17:23]) {
        
        svg(paste0('plots/others/index-',i,'-linear.svg'), width = 8.57, height = 5.79)
        
        df <- index_all[,c("celltype","seurat_clusters",`i`)]
        colnames(df)[4] <- "Protein"
        pl <-  ggplot(df, aes(factor(celltype, levels=levels(metadata$celltype)), Protein, fill=celltype)) +
                geom_violin(scale = "width", lwd=0.5, draw_quantiles = 0.5) +
                #coord_flip() +
                #theme_pubr() +
                theme_pubclean() +
                labs(y="Expression (A.U.)",x=element_blank(), title = i) +
                scale_fill_manual(values = colors_many, guide=F) + 
                #scale_y_continuous(limits=c(0,9), breaks = c(0,5)) +
                theme(
                        strip.text.x = element_text(
                                size = 12, color = "black", face="bold"
                        ))
        print(pl)
        dev.off()
        
}


