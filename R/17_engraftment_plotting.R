library(readxl)
library(tidyverse)
library(tidyquant)
library(ggpubr)
library(ggbeeswarm)

#load data
df <- read_excel(file.path("data","engraftment","table_S1.xlsx")) %>% 
  as.data.frame()
rownames(df) <- df[,1]
df <- df[,-1]
colnames(df) <- paste0("Pat_", colnames(df))

#subset for the relevant patients
df2 <- df[c("LM engraftment IBA1 [%]","LM engraftment CD206 [%]","LM engraftment SIGLEC1 [%]","Pv engraftment IBA1 [%]", "Pv engraftment CD206 [%]", "Pv engraftment SIGLEC1 [%]"),c("Pat_1", "Pat_12","Pat_14")] %>% 
  rownames_to_column(var="Parameter") %>% 
  mutate(Antigen=factor(rep(c("IBA1","CD206","SIGLEC1"),2),levels = c("IBA1","CD206","SIGLEC1"))) %>% 
  pivot_longer(-c(Parameter, Antigen), names_to = "Pat_ID" , values_to = "percent") %>% 
  mutate(Compartment=case_when(
    grepl("LM",.$Parameter) ~ "Leptomeninges",
    T ~ "Perivascular space"
  ),
  percent=as.numeric(percent))

#dot plot
my_comparisons <- list( c("IBA1", "CD206"), c("CD206", "SIGLEC1"), c("IBA1", "SIGLEC1") )

df2 %>% 
  ggplot(aes(Antigen, percent, fill=Antigen)) +
  geom_jitter(shape=21, size=12, width = .2, height = 0) +
  stat_summary(fun=mean, geom="crossbar", width=0.7) +
  facet_wrap(~ Compartment, ncol = 1) +
  expand_limits(y=c(0,100)) +
  theme_linedraw() +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test",
                     size=7) +
  labs(y="% of double positive cells") +
  theme(panel.grid = element_blank(),
        legend.position = "None",
        text = element_text(size = 25)) +
  scale_fill_brewer(palette = "Pastel2")
  
ggsave(file.path("plots", "others", "engraftment","engraftment_dot_plot.pdf"), height = 10, width=5, useDingbats=F)
