library(readxl)
library(tidyverse)
library(rstatix)
library(DescTools)
library(ggplot2)
library(multcomp)
library(gridExtra)

exons <- readRDS(file="PTCH1_exons_rel.rds")

ggplot(exons, aes(x=Day, y=-ratio, color=Line)) + 
  geom_point() + 
  geom_smooth(method="lm", alpha=0.2) +
  theme_classic() +
  scale_y_continuous(expand=expansion(mult=0.2))+
  scale_color_manual(values = c("#AED6F1", "#FFC300"))

# Drug targets
#packages:
library(readxl)
library(tximeta)
library(SummarizedExperiment)
library(org.Dm.eg.db)
library(DESeq2)
library(tidyverse)
library("vsn")
library(pheatmap)
library("RColorBrewer")
library("genefilter")
library(fgsea)
library(org.Hs.eg.db)
library(tibble)
library(dplyr)
library(ggplot2)
library(magrittr)
library("genefilter")
library(gridExtra)
library(grid)
library(clusterProfiler)
library(EnsDb.Hsapiens.v86)
library("AnnotationHub")
library("AnnotationDbi")
library(ensembldb)
library(ggvenn)
library("vsn")
library(GSVA)
library(msigdbr)
library(limma)
library(ggrepel)
library(cowplot)
library(clusterProfiler)
library(enrichplot)

setwd("/Users/max/Documents/PTCH1_manuscript/data_code/source")

# Figure 3A
dds3 <- readRDS("dds3.rds")
vsd3 <- readRDS("vsd3.rds")

FDA_drugs <- c("ERBB4",  "YAP1", "HDAC1", "HDAC9", "HMMR", "TOP2A",  "IDH1", "SMO", "SQLE", "HMGCR", "PCSK9", "FDFT1")
rownames(vsd3) <- mcols(dds3)$gene_name
pheatmap(assay(vsd3)[FDA_drugs, ], scale="row", cellwidth=12, cellheight=12, treeheight_col=4, treeheight_row=4, border_color="white", show_colnames = F, legend=F)

