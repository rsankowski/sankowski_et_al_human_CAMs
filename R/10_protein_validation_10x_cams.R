library(tidyverse)
library(readxl)
library(broom)
library(ggbeeswarm)
library(ggpubr)
library(RColorBrewer)

source("R/functions.R")

load(file.path("data","protein_validation_cams","protein_valid_cumulative_counts.RData")) 

prot_all$Compartment2 <- ifelse(prot_all$Localization == "CP Macs", "cpMΦ", 
                                ifelse(prot_all$Localization == "Dura", "dMΦ", 
                                       ifelse(prot_all$Localization == "CP Epi", "CP_epi",
                                              ifelse(prot_all$Compartment == "PV", "pvMΦ",
                                                     ifelse(prot_all$Localization == "Men", "mMΦ",prot_all$Compartment))))) %>%
  factor(levels = c("Micr", "pvMΦ", "cpMΦ", "CP_epi", "mMΦ", "dMΦ"))

prot_all$Antigen <- factor(prot_all$Antigen, levels = c("CD206_IBA1", "SIGLEC1_CD206", "SPP1_Iba1","S100A6_CD206","CD1C_IBA1","SPP1_Iba1"))

#define comparisons
my_comparisons <- apply(expand.grid("pvMΦ",c("Micr",  "cpMΦ", "CP_epi", "mMΦ", "dMΦ"), stringsAsFactors = F), 1, list)
my_comparisons <- lapply(my_comparisons, unlist)
my_comparisons <- lapply(my_comparisons, unname)

#define colors
cols <- c(brewer.pal(4,"Set1"))
my_cols <- c()
my_cols[5:6] <- cols[3:2]
my_cols[1:2] <- colors_fig[c(1,3)]
my_cols[3:4] <- c("#CB4154","#EE7942")


walk(unique(prot_all$Antigen), function(x) {
  plt <- prot_all %>%
  filter(Antigen == x) %>%
  ggplot(aes(x= Compartment2, y=median_ratio*100, fill=Compartment2)) +
  geom_quasirandom(size=15, pch=21) +
  stat_summary(fun = "median", geom="crossbar", color="black") +
  expand_limits(y=c(0,1)) +
  theme_linedraw() +
  theme(legend.position = NaN,
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size=30)) +
  scale_fill_manual(values = my_cols) +
  stat_compare_means(comparisons = my_comparisons, hide.ns = T, size=10) +
  labs(y="% pos. cells", x=element_blank()) +
    expand_limits(y=c(0,100))

  print(plt)
  
ggsave(file.path("plots","others","10x_cams",paste("dotplot_ctrl_iba1",x,".pdf", sep="_")), useDingbats=F, height = 8, width=8)
})

