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

# PCA control + hets:
pca_ua_r2adb <- prcomp(t(assay(vsd3)))

pca_ua_r2adbdf <- data.frame(pca_ua_r2adb$x)
rownames(pca_ua_r2adbdf) <- colnames(assay(vsd3))
pca_ua_r2adbdf <- cbind(pca_ua_r2adbdf, vsd3$genotype, vsd3$genotype_corrected)
colnames(pca_ua_r2adbdf)[13] <- c("genotype")
pcaplot_adj <- ggplot(pca_ua_r2adbdf, aes(x=PC1, y=PC2, colour=genotype)) + 
  geom_point(size=3) +
  theme_classic()+
  labs(x="PC1 (")
ggtitle("PCA")
pcaplot_adj

## Heatmap
setwd("/Users/max/Documents/PTCH1_manuscript/data_code/source")

rownames(vsd3) <- mcols(dds3)$gene_name

##
RL_genes <- c("PAX6", "ZIC1", "ZIC2", "PCNA", "MKI67", "DLGAP5", "CCND1", 
              "NEUROD1", "NEUROD2", "CNTN1", "CNTN2", "UNC5C", 
              "EOMES", "CALB2", "GRM1","TBR1", "GLI1", "GLI2", "GLI3", "PTCH1")

ca <- data.frame(vsd3$genotype)
rownames(ca) <- colnames(vsd3)
colnames(ca) <- "Genotype"

pheatmap(assay(vsd3)[RL_genes, ], scale="row",
         annotation_col = ca, cluster_rows = F, cluster_cols = F,  show_colnames = F, cellwidth = 10, cellheight = 10,
         treeheight_col=4, treeheight_row=4, border_color="white", legend_labels = c("", ""))

# Euclidean distances
setwd("/Users/max/Documents/PTCH1_manuscript/data_code/source")
gene_names <- read_RDS(file="gene_names.rds")
vsd4 <- readRDS("vsd4.rds")
aldinger_genes <- read_csv2("Aldinger_genes.csv")

vsd4_2 <- vsd4
rownames(vsd4_2) <- gene_names
sampleDists <- dist(t(assay(vsd4_2)[rownames(vsd4_2)%in%aldinger_genes$Gene, ]))
eucl <- as.matrix(sampleDists)[51:61, 1:50]
eanno <- data.frame(colData(vsd4_2))

eanno <- eanno[1:50, c(4, 5)]
eanno$Region <- factor(eanno$Region, levels=c("RL", "EGL", "PCL", "Bulk"), ordered=T)
eanno <- eanno[with(eanno, order(Region, age_pcw)), ]
eucl <- eucl[ , rownames(eanno)]
eanno2 <- data.frame(colData(vsd4))[c(51:53, 56:58, 54, 55, 59, 60, 61), c(4, 5)]
cannohu <- rbind(eanno, eanno2)
colors_sd2 <- colorRampPalette(brewer.pal(10, "RdYlBu"))(11)
sampdistsph <- pheatmap(eucl, cluster_cols = F, scale="row", annotation_col=eanno, border_color = "white",
                        color = colors_sd2, show_rownames = T, show_colnames=T, cellwidth = 10, cellheight = 10, legend=F, annotation_legend = TRUE)

## Dot/treeplot combination
## Dot/treeplot combination
#Finding unique genes:
##
het_res <- lfcShrink(dds3, coef = "genotype_Het_vs_WT")
het_res$gene <- mcols(dds3)$gene_name

hom_res <- lfcShrink(dds2, coef = "genotype_Hom_vs_WT")
hom_res$gene <- mcols(dds2)$gene_name

het_res_sig <- subset(het_res, het_res$padj<0.05)
hom_res_sig <- subset(hom_res, hom_res$padj<0.05)

#Finding unique genes:
outersect1 <- function(x, y) {setdiff(x, y)}
##
hetsigdf <- data.frame(het_res_sig)
hetsigdf$type <- rep("ht", nrow(hetsigdf))
hetsigdf$rn <- rownames(hetsigdf)
homsigdf <- data.frame(hom_res_sig)
homsigdf$type <- rep("hm", nrow(homsigdf))
homsigdf$rn <- rownames(homsigdf)
bothsigdf <- merge(hetsigdf, homsigdf, by="rn", all=TRUE)
bothsigdf$oposite <- ifelse(bothsigdf$log2FoldChange.x>0&bothsigdf$log2FoldChange.y<0, "a",
                            ifelse(bothsigdf$log2FoldChange.x<0&bothsigdf$log2FoldChange.y>0, "a","b"))

oposites <- subset(bothsigdf, bothsigdf$oposite=="a")

ht_unique <- c(ht_only, oposites$rn)
ht_unique_res <- subset(hetsigdf, rownames(hetsigdf) %in% ht_unique)
ht_unique_up <- subset(ht_unique_res, ht_unique_res$log2FoldChange>0)

universe <- readRDS("universe.rds")

# Construct differentially upregulated pathway matrix
hu_bp <- clusterProfiler::enrichGO(gene= rownames(ht_unique_up),
                                   universe      = universe,
                                   OrgDb         = "org.Hs.eg.db",
                                   keyType = "ENSEMBL",
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)

#dotplot(hu_bp, showCategory=30)
setReadable(hu_bp, 'org.Hs.eg.db', 'ENTREZID') %>%
  pairwise_termsim(.) %>% treeplot(., showCategory = 30, nCluster=5)
#Generated plots above were combined together in Adobe illustrator

## Volcanoplot
genes2 <- c("PAX6", "LHX5", "ANOS1", "EPHA4", "GDF10", "GBX2", "CEP290", "TMEM67", "WLS", "WNT7A")
data.frame(ht_unique_res)%>%
  subset(gene%in%genes2) ->imps

data.frame(ht_unique_res)%>%
  subset(log2FoldChange>0) %>% 
  ggplot(., aes(x = log2FoldChange,
                y = -log10(padj)))+
  geom_point(alpha=0.3, aes(size=baseMean), color="grey") +
  geom_point(alpha=0.5, data=imps, aes(size=baseMean), 
             color="orange")+
  geom_label_repel(data = imps,   
                   aes(label = gene),
                   max.overlaps = Inf) +
  theme_classic()
# Figure colors were further amended in adobe illustrator
