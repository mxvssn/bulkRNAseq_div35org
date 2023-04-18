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

setwd("/Users/max/Documents/PTCH1_manuscript/data_code/source")

# Figure 3A
dds1 <- readRDS("dds1.rds")
vsd1 <- readRDS("vsd1.rds")

geno_cols <- c("#AED6F1", "#FFC300", "#FF5733")

topVarGenes5000 <- head(order(rowVars(assay(vsd1)), decreasing = TRUE), 1000)
r1_var <-vsd1[topVarGenes5000, ]
pca1 <- prcomp(t(assay(r1_var)))

pca1df <- data.frame(pca1$x)
rownames(pca1df) <- colnames(assay(vsd1))
pca1df <- cbind(pca1df, vsd1$genotype)
colnames(pca1df)[9] <- c("genotype")
pcaplot1 <- ggplot(pca1df, aes(x=PC1, y=PC2, colour=genotype, label=rownames(pca1df))) + 
  geom_point(size=3, alpha=0.5) +
  theme_classic()+
  scale_color_manual(values=geno_cols)+
  labs(x="PC1 (74% of variance explained)", y="PC2 (12% of variance explained)", col="Genotype")
pcaplot1

# Differential gene expression testing control vs PTCH1 homozygous
dds2 <- readRDS("dds2.rds")
vsd2 <- readRDS("vsd2.rds")

hom_res <- lfcShrink(dds2, coef = "genotype_Hom_vs_WT")
hom_res$gene <- mcols(dds2)$gene_name

summary(hom_res$baseMean>50&hom_res$padj<0.05, na.rm=TRUE)

hom_res_sig <- subset(hom_res, hom_res$baseMean>50&hom_res$padj<0.05)

#Figure 3B
#Volcanoplot homozygotes
top_hom <- hom_res_sig[order(hom_res_sig$padj), ] %>% .[1:500, ]
top_hom_up <- top_hom[order(top_hom$log2FoldChange), ] %>% .[1:10, ]
top_hom_down <- top_hom[order(top_hom$log2FoldChange, decreasing=T), ] %>% .[1:10, ]
top_hom2 <- rbind(top_hom_up, top_hom_down)

hom_res_mat <- data.frame(hom_res)

hom_lab <- subset(hom_res_mat, hom_res_mat$gene %in% top_hom2$gene)

cols_hom <- c("up" = "#AED6F1", "down" = "#FF5733", "ns" = "grey") 

hom_res_mat <- data.frame(hom_res_sig) %>%
  mutate(gene_type= case_when(log2FoldChange>=0.5849625 & padj <=0.05 ~ "Up", 
                              log2FoldChange<= c(-0.5849625) & padj<=0.05 ~ "Down", 
                              TRUE ~ "ns"))


hom_volc <- ggplot(data=hom_res_mat, aes(x = log2FoldChange,
                                         y = -log10(padj)))+
  geom_point(alpha=0.3, aes(color=gene_type)) +
  geom_point(data = hom_lab,
             size = 2,
             shape = 21,
             fill = "firebrick",
             colour = "black")  +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             alpha=0.5) +
  geom_vline(xintercept = c(-0.585, 0.585),
             linetype = "dashed",
             alpha=0.5)+
  geom_label_repel(data = hom_lab,   
                   aes(label = gene),
                   max.overlaps = Inf) +
  scale_color_manual(values = c("#AED6F1", "grey", "#FF5733"))+
  labs(x = "log2(fold change)",
       y = "-log10(BH adjusted P-value)",
       colour = "Expression \nchange") +
  theme_classic()
hom_volc

# Gene expression region genes
mat_hom <- assay(vsd2)
rownames(mat_hom) <- mcols(dds2)$gene_name

#Specify genes
FP <- c("SHH","NTN1", "FOXA2", "NKX6-1", "NKX2-2")
RP <- c("POU4F1", "GDF7", "RSPO1", "PAX6", "WNT3A")
FB <- c("FOXG1", "NKX2-1", "SIX3")
HB <- c("GBX2", "LMX1A", "ATOH1")
fig3_genes <- c(FP, RP, FB, HB)

#Row annotation
ra_fig3 <- data.frame(
  factor(
    c(rep("Ventral", 5), 
      rep("Dorsal", 5), 
      rep("Forebrain", 3), 
      rep("Hindbrain", 3)), 
    levels=c("Ventral", "Dorsal", "Forebrain", "Hindbrain"), ordered=T))

colnames(ra_fig3)[1] <- "Region"
rownames(ra_fig3) <- fig3_genes

#Column annotation
ca2 <- data.frame(vsd2$genotype)
rownames(ca2) <- colnames(vsd2)
colnames(ca2)[1] <- "Genotype"

#Heatmap
pheatmap(mat_hom[fig3_genes,], scale="row", main="Brain regions", cluster_rows=F, show_rownames=T, show_colnames = F,
         cluster_cols=FALSE, annotation_row=ra_fig3, annotation_col = ca2, gaps_row=c(10,16),
         cellwidth = 12, cellheight = 12, border_color = "white", legend=F)