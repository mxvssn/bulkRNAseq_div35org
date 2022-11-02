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

#set wd
setwd("/Volumes/Partitie4/RNAseq_analysis_2/")
load("/Volumes/Partitie4/RNAseq_analysis_2/RNAseqR2_5.RData")


#Load dataset
test <- read_xlsx("/Volumes/Partitie4/RNAseq_analysis_2/samples4.xlsx")
#file.exists(test$files)
#View(test)
#Transcript quantification import and generation of RangedSummarizedExperiment
# Ensembl - Homo sapiens - release 97 is used as transcriptome
se <- tximeta(test)
# assayNames(se)
# rowRanges(se)
# seqinfo(se)
# edb <- retrieveDb(se)
# class(edb)

# Generate count matrix for gene counts
gse <- summarizeToGene(se)

#Generate dataset which including only WT and Het:
gse_r2ad <- gse[, -c(7, 13:23)]
dds_r2ad <- DESeqDataSet(gse_r2ad, ~ batch + genotype)
dds_r2ad$genotype <- relevel(dds_r2ad$genotype, ref="WT")

# data exploration ----
dds_r2ad <- DESeq(dds_r2ad,
                 fitType = "parametric",
                 minReplicatesForReplace = Inf)

#determine which fit is most approriate:
par_r2ad <- estimateDispersions(dds_r2ad, fitType = "parametric")
residualspar_r2ad <- abs(log(mcols(par_r2ad)$dispGeneEst) - log(mcols(par_r2ad)$dispFit))
loc_r2ad <- estimateDispersions(dds_r2ad, fitType = "local")
residualsloc_r2ad <- abs(log(mcols(loc_r2ad)$dispGeneEst) - log(mcols(loc_r2ad)$dispFit))
summary(residualspar_r2ad) # median is 1.013
summary(residualsloc_r2ad) # median is 0.854
# Maintain parametric fit!

saveRDS(dds_r2ad, file="dds_r2ad.rds")

#Varient stabalizing transformation
vsd_r2ad <- vst(dds_r2ad, blind=FALSE, fitType = "parametric")
saveRDS(vsd_r2ad, file="vsd_r2ad.rds")

# Different normalisation methods
rld <- rlog(dds_r2ad, blind=FALSE)
ntd <- normTransform(dds_r2ad)
SDnorm <- meanSdPlot(assay(ntd))
SDrld <- meanSdPlot(assay(rld))
SDvsd <- meanSdPlot(assay(vsd_r2ad))
SDnorm2 <- SDnorm$gg + ggtitle("mean SD normalized log2(n+1)") + scale_fill_gradient(low = "blue", high = "darkred") + scale_y_continuous(limits = c(0, 7))
leg <- get_legend(SDnorm2)
SDnorm2 <- SDnorm$gg + ggtitle("mean SD normalized log2(n+1)") + scale_fill_gradient(low = "blue", high = "darkred") + scale_y_continuous(limits = c(0, 7))  + theme_classic() +
  theme(legend.position="none")
SDrld2 <- SDrld$gg + ggtitle("mean SD normalized rlog") + scale_fill_gradient(low = "blue", high = "darkred") + scale_y_continuous(limits = c(0, 7))  + theme_classic() +
  theme(legend.position="none")
SDvsd2 <- SDvsd$gg + ggtitle("mean SD normalized vst") + scale_fill_gradient(low = "blue", high = "darkred") + scale_y_continuous(limits = c(0, 7))  + theme_classic() +
  theme(legend.position="none")
plot_grid(SDnorm2, SDrld2, SDvsd2, leg, ncol=4, rel_widths = c(3/10, 3/10, 3/10, 1/10))

SDnorm$gg + ggtitle("mean SD normalized log2(n+1)") + scale_fill_gradient(low = "blue", high = "darkred") + scale_y_continuous(limits = c(0, 7)) + theme_classic() +
  theme(legend.position="none")

par(mfrow=c(1,2), mar=c(4,4,2,1))
countsraw <- boxplot(log2(counts(dds_r2ad, normalized=FALSE)+1), las=2, main="Raw counts", col=c(rep("lightblue", 7), rep("orange", 5), rep("red", 3)))
countsntd <- boxplot(assay(ntd),las=2,  main="Normalized counts: log2(n+1)", col=c(rep("lightblue", 4), rep("orange", 2), rep("red", 4)))
countsrld <- boxplot(assay(rld),las=2,  main="Normalized counts: rlog", col=c(rep("lightblue", 4), rep("orange", 2), rep("red", 4)))
countsvsd_r2ad <- boxplot(assay(vsd_r2ad),las=2,  main="Normalized counts: vst", col=c(rep("lightblue", 7), rep("orange", 5), rep("red", 2)))
par(mfrow=c(1,1), mar=c(4,4,2,1))

# PCA:
pca_ua_r2ad <- prcomp(t(assay(vsd_r2ad)))

dims <- data.frame(PC= factor(paste0("PC",1:15), levels=c(paste0("PC",1:15)),ordered=T),
                   var_explained=(pca_ua_r2ad$sdev)^2/sum((pca_ua_r2ad$sdev)^2))
head(dims, 4)
elbow<- dims[1:15,]
screeplot <- ggplot(elbow, aes(x=PC,y=var_explained, group=1))+
  geom_point(size=4)+
  geom_line()+
  labs(title="Scree plot, wt and ht only, unadjusted") +
  ylab("Variance explained")+
  theme_classic()
screeplot


pca_ua_r2addf <- data.frame(pca_ua_r2ad$x)
rownames(pca_ua_r2addf) <- colnames(assay(vsd_r2ad))
pca_ua_r2addf <- cbind(pca_ua_r2addf, vsd_r2ad$genotype, vsd_r2ad$genotype_corrected)
colnames(pca_ua_r2addf)[13] <- c("genotype")
pcaplot<- ggplot(pca_ua_r2addf, aes(x=PC1, y=PC2, colour=genotype, label=rownames(pca_ua_r2addf))) + 
  geom_point(size=3) +
  geom_text_repel(max.overlaps = Inf) +
  theme_classic()+
  ggtitle("PCA, wt and ht only, unadjusted")
pcaplot

#Batch correction:
vsd_r2adb <- vsd_r2ad
assay(vsd_r2adb) <- limma::removeBatchEffect(assay(vsd_r2adb), vsd_r2adb$batch)

# PCA:
pca_ua_r2adb <- prcomp(t(assay(vsd_r2adb)))

dims <- data.frame(PC= factor(paste0("PC",1:11), levels=c(paste0("PC",1:11)),ordered=T),
                   var_explained=(pca_ua_r2adb$sdev)^2/sum((pca_ua_r2adb$sdev)^2))
head(dims, 4)
elbow<- dims[1:15,]
screeplot_adj <- ggplot(elbow, aes(x=PC,y=var_explained, group=1))+
  geom_point(size=4)+
  geom_line()+
  labs(title="Scree plot, wt and ht only, limma adjusted") +
  ylab("Variance explained")+
  theme_classic()
screeplot_adj


pca_ua_r2adbdf <- data.frame(pca_ua_r2adb$x)
rownames(pca_ua_r2adbdf) <- colnames(assay(vsd_r2adb))
pca_ua_r2adbdf <- cbind(pca_ua_r2adbdf, vsd_r2adb$genotype, vsd_r2adb$genotype_corrected)
colnames(pca_ua_r2adbdf)[13] <- c("genotype")
pcaplot_adj <- ggplot(pca_ua_r2adbdf, aes(x=PC1, y=PC2, colour=genotype)) + 
  geom_point(size=3) +
  theme_classic()+
  labs(x="PC1 (")
  ggtitle("PCA")
pcaplot_adj

# Assessing sample similarity
sampleDists <- dist(t(assay(vsd_r2adb)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd_r2adb$genotype, vsd_r2adb$line, vsd_r2adb$batch, sep="-")
colnames(sampleDistMatrix) <- NULL
sampleDistplot <- pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)


## DEG for WT vs Het
resultsNames(dds_r2ad)

dds_r2ad <- readRDS("/Volumes/Partitie4/RNAseq_analysis_2dds_r2ad.rds")

res_ht_r2ad <- results(dds_r2ad, name = "genotype_Het_vs_WT", alpha=0.05)
summary(res_ht_r2ad$padj<0.05, na.rm=TRUE) #3821 genes
shr_ht_r2ad <- lfcShrink(dds_r2ad, coef = "genotype_Het_vs_WT")
summary(shr_ht_r2ad$padj<0.05& shr_ht_r2ad$baseMean>50, na.rm=TRUE) #3821 genes
shr_ht_r2ad$gene <- mcols(dds_r2ad)$gene_name
write.csv(shr_ht_r2ad,  file="shr_ht_r2ad.csv")

res_ht_r2ad_sig <- subset(res_ht_r2ad, res_ht_r2ad$padj<0.05)
shr_ht_r2ad_sig <- subset(shr_ht_r2ad, shr_ht_r2ad$padj<0.05& shr_ht_r2ad$baseMean>50)
summary(shr_ht_r2ad_sig$log2FoldChange>0.5849625, na.rm=T)

# PCA with DGEs:
vsd_ht_dge <- subset(vsd_r2adb, rownames(vsd_r2adb)%in%rownames(shr_ht_r2ad_sig))
head(order(rowVars(assay(CBTNvsd)), decreasing = TRUE), 600)
topVarGenes5000 <- head(order(rowVars(assay(vsd_r2adb)), decreasing = TRUE), 5000)
vsd_ht_var <-vsd_r2adb[topVarGenes5000, ]
pca_ua_ht_var <- prcomp(t(assay(vsd_ht_var)))

dims <- data.frame(PC= factor(paste0("PC",1:11), levels=c(paste0("PC",1:11)),ordered=T),
                   var_explained=(pca_ua_ht_var$sdev)^2/sum((pca_ua_ht_var$sdev)^2))
head(dims, 4)
elbow<- dims[1:15,]
screeplot_adj <- ggplot(elbow, aes(x=PC,y=var_explained, group=1))+
  geom_point(size=4)+
  geom_line()+
  labs(title="Scree plot, wt and ht only, limma adjusted") +
  ylab("Variance explained")+
  theme_classic()
screeplot_adj


pca_ua_ht_vardf <- data.frame(pca_ua_ht_var$x)
rownames(pca_ua_ht_vardf) <- colnames(assay(vsd_ht_var))
pca_ua_ht_vardf <- cbind(pca_ua_ht_vardf, vsd_ht_var$genotype, vsd_ht_var$batch)
colnames(pca_ua_ht_vardf)[c(12, 13)] <- c("genotype", "batch")
pcaplot_adj <- ggplot(pca_ua_ht_vardf, aes(x=PC1, y=PC2, colour=genotype)) + 
  geom_point(size=3, aes(shape=factor(batch))) +
  theme_classic()+
  labs(x="PC1 (52.7% eplained)", y="PC2 (12.2% explained)", title="PCA", subtitle="Top 5000 most variable genes")
pcaplot_adj

#MA plots
DESeq2::plotMA(res_ht_r2ad, ylim=c(-3,3), cex=.8, main="WT vs Het, p-value < 0.05")
DESeq2::plotMA(shr_ht_r2ad, ylim=c(-3,3), cex=.8, main="WT vs Het, lfcShrink, p-value < 0.05")

## DEGs:
ca <- data.frame(dds_r2ad$genotype)
ca$batch <- dds_r2ad$batch
rownames(ca) <- colnames(dds_r2ad)
colnames(ca) <- c("Genotype", "Batch")
cano <- c(rep("WT", 6), rep("Het", 5))
pheatmap(assay(vsd_r2adb)[rownames(shr_ht_r2ad_sig), ], scale="row", show_rownames = F, annotation_col=ca, show_colnames = F,  treeheight_col = 5, treeheight_row = 5)

### Now homs only:
#Generate dataset which including only WT and Het:
gse_r2ahm <- gse[, -c(4:12, 16:23)]
dds_r2ahm <- DESeqDataSet(gse_r2ahm, ~  genotype)
dds_r2ahm$genotype <- relevel(dds_r2ahm$genotype, ref="WT")

# data exploration ----
dds_r2ahm <- DESeq(dds_r2ahm,
                  fitType = "parametric",
                  minReplicatesForReplace = Inf)

# saveRDS(dds_r2ahm, file="dds_r2ahm.rds")
dds_r2ahm <- readRDS("dds_r2ahm.rds")

#Varient stabalizing transformation
vsd_r2ahm <- vst(dds_r2ahm, blind=FALSE, fitType = "parametric")
# saveRDS(vsd_r2ahm, file="vsd_r2ahm.rds")

res_hm_r2ahm <- results(dds_r2ahm, name = "genotype_Hom_vs_WT", alpha=0.05)
summary(res_hm_r2ahm$padj<0.05, na.rm=TRUE) #4608 genes
shr_hm_r2ahm <- lfcShrink(dds_r2ahm, coef = "genotype_Hom_vs_WT")
summary(shr_hm_r2ahm$padj<0.05, na.rm=TRUE) #4579 genes
summary(shr_hm_r2ahm$baseMean>50&shr_hm_r2ahm$padj<0.05, na.rm=TRUE)
shr_hm_r2ahm$gene <- mcols(dds_r2ahm)$gene_name
write.csv(shr_hm_r2ahm,  file="shr_hm_r2ahm.csv")


res_hm_r2_sig <- subset(res_hm_r2ahm, res_hm_r2ahm$padj<0.05)
shr_hm_r2_sig <- subset(shr_hm_r2ahm, shr_hm_r2ahm$baseMean>50&shr_hm_r2ahm$padj<0.05)


## DEGs homs:
cah <- data.frame(dds_r2ahm$genotype)
rownames(cah) <- colnames(dds_r2ahm)
colnames(cah) <- c("Genotype")
canoh <- c(rep("WT", 3), rep("Hom", 3))

top_hom <- shr_hm_r2_sig[order(shr_hm_r2_sig$padj), ] %>% .[1:500, ]
top_hom_up <- top_hom[order(top_hom$log2FoldChange), ] %>% .[1:10, ]
top_hom_down <- top_hom[order(top_hom$log2FoldChange, decreasing=T), ] %>% .[1:10, ]
top_hom2 <- rbind(top_hom_up, top_hom_down)
pheatmap(assay(vsd_r2ahm)[rownames(top_hom2), ], scale="row", show_rownames = F,
         annotation_col = cah, cellwidth = 20)


#Finding unique genes:
outersect1 <- function(x, y) {setdiff(x, y)}
##
hetsigdf <- data.frame(shr_ht_r2ad_sig)
hetsigdf$type <- rep("ht", nrow(hetsigdf))
hetsigdf$rn <- rownames(hetsigdf)
homsigdf <- data.frame(shr_hm_r2_sig)
homsigdf$type <- rep("hm", nrow(homsigdf))
homsigdf$rn <- rownames(homsigdf)
bothsigdf <- merge(hetsigdf, homsigdf, by="rn", all=TRUE)
bothsigdf$oposite <- ifelse(bothsigdf$log2FoldChange.x>0&bothsigdf$log2FoldChange.y<0, "a",
                            ifelse(bothsigdf$log2FoldChange.x<0&bothsigdf$log2FoldChange.y>0, "a","b"))

oposites <- subset(bothsigdf, bothsigdf$oposite=="a")

ht_only <- outersect1(rownames(shr_ht_r2ad_sig), rownames(shr_hm_r2_sig))
hetsigdf <- data.frame(shr_ht_r2ad_sig)
hetsigdf$type <- rep("ht", nrow(hetsigdf))
hetsigdf$rn <- rownames(hetsigdf)
homsigdf <- data.frame(shr_hm_r2_sig)
homsigdf$type <- rep("hm", nrow(homsigdf))
homsigdf$rn <- rownames(homsigdf)
bothsigdf <- merge(hetsigdf, homsigdf, by="rn", all=TRUE)
bothsigdf$oposite <- ifelse(bothsigdf$log2FoldChange.x>0&bothsigdf$log2FoldChange.y<0, "a",
                            ifelse(bothsigdf$log2FoldChange.x<0&bothsigdf$log2FoldChange.y>0, "a","b"))

oposites <- subset(bothsigdf, bothsigdf$oposite=="a")
ht_unique <- c(ht_only, oposites$rn)
length(ht_unique) # -> 1577 genes specifically differentially expressed in het comp wt
ht_unique_res2 <- shr_ht_r2ad[ht_unique, ]

pheatmap(assay(vsd_r2adb)[ht_unique, ], scale="row", show_rownames = F,
         main="DEGs uniquely upregulated in heterozygous", annotation_col = ca,
         cellwidth = 20)

ht_unique_res <- subset(shr_ht_r2ad, rownames(shr_ht_r2ad) %in% ht_unique)
# View(data.frame(ht_unique_res_hbm))
# ggplot(data.frame(ht_unique_res), aes(y=baseMean))+geom_histogram(bins=500)+scale_y_log10()
# ht_unique_res_hbm <- subset(ht_unique_res, ht_unique_res$baseMean>1000)

ht_unique_genes <- rownames(ht_unique_res)
View(data.frame(ht_unique_genes))

write.csv(ht_unique_genes,  file="ht_unique_genes.csv")

save.image(file = "RNAseqR2_5.RData")

### Which GO's are enriched:

ht_unique_up <- subset(ht_unique_res, ht_unique_res$log2FoldChange>0)
ht_unique_down <- subset(ht_unique_res, ht_unique_res$log2FoldChange<0)
universe <- rownames(shr_ht_r2ad)

GOup <- clusterProfiler::enrichGO(gene= rownames(shr_ht_r2ad_sig[shr_ht_r2ad_sig$log2FoldChange>0, ]),
                                   universe      = universe,
                                   OrgDb         = "org.Hs.eg.db",
                                   keyType = "ENSEMBL",
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)

dotplot(GOsig, showCategory=10) + ggtitle("dotplot for GO enirchement")


go_up <- clusterProfiler::enrichGO(gene= ht_unique_up,
                universe      = universe,
                OrgDb         = "org.Hs.eg.db",
                keyType = "ENSEMBL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

go_down <- clusterProfiler::enrichGO(gene= ht_unique_down,
                                   universe      = universe,
                                   OrgDb         = "org.Hs.eg.db",
                                   keyType = "ENSEMBL",
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)


go_up_cc <- clusterProfiler::enrichGO(gene= ht_unique_up,
                                   universe      = universe,
                                   OrgDb         = "org.Hs.eg.db",
                                   keyType = "ENSEMBL",
                                   ont           = "CC",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)

go_down_cc <- clusterProfiler::enrichGO(gene= ht_unique_down,
                                     universe      = universe,
                                     OrgDb         = "org.Hs.eg.db",
                                     keyType = "ENSEMBL",
                                     ont           = "CC",
                                     pAdjustMethod = "BH",
                                     pvalueCutoff  = 0.01,
                                     qvalueCutoff  = 0.05,
                                     readable      = TRUE)

go_up_mf <- clusterProfiler::enrichGO(gene= ht_unique_up,
                                      universe      = universe,
                                      OrgDb         = "org.Hs.eg.db",
                                      keyType = "ENSEMBL",
                                      ont           = "MF",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff  = 0.01,
                                      qvalueCutoff  = 0.05,
                                      readable      = TRUE)

go_down_mf <- clusterProfiler::enrichGO(gene= ht_unique_down,
                                        universe      = universe,
                                        OrgDb         = "org.Hs.eg.db",
                                        keyType = "ENSEMBL",
                                        ont           = "MF",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff  = 0.01,
                                        qvalueCutoff  = 0.05,
                                      readable      = TRUE)

View(data.frame(go_up))
View(data.frame(go_down))
View(data.frame(go_up_cc))
View(data.frame(go_down_cc))
View(data.frame(go_up_mf))
View(data.frame(go_down_mf))

library(DOSE)
data(geneList)
View(geneList)


ht_unique_up2 <- rownames(ht_unique_res[ht_unique_res$baseMean>50&ht_unique_res$log2FoldChange>0.5, ])
ht_unique_down2 <- rownames(ht_unique_res[ht_unique_res$baseMean>50&ht_unique_res$log2FoldChange<c(-0.5), ])
length(ht_unique_down2)

go_up <- clusterProfiler::enrichGO(gene= ht_unique_up2,
                                   universe      = universe,
                                   OrgDb         = "org.Hs.eg.db",
                                   keyType = "ENSEMBL",
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
nrow(go_up)

go_down <- clusterProfiler::enrichGO(gene= ht_unique_down2,
                                   universe      = universe,
                                   OrgDb         = "org.Hs.eg.db",
                                   keyType = "ENSEMBL",
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
nrow(go_down)
p1 <- dotplot(go_up) + scale_color_viridis_c() + theme_classic() + 
  labs(title="Enriched gene sets enriched of uniquely upregulated genes in heterozygous", x="Gene ratio", y="")

p2 <- dotplot(go_down) + scale_color_viridis_c() + theme_classic() +
  labs(title="Enriched gene sets enriched of uniquely downregulated genes in heterozygous", x="Gene ratio", y="")

grid.arrange(p1, p2, ncol=1)

ggplot(go_up[order(go_up$GeneRatio), ], aes(x=GeneRatio, y=Description, size=Count, color=p.adjust))+
  geom_point()

p3 <- setReadable(go_up, 'org.Hs.eg.db', 'ENTREZID') %>% enrichplot::pairwise_termsim(.) %>% treeplot(.) + scale_color_viridis_c() + scale_fill_brewer(palette="Set1")
p4 <- setReadable(go_down, 'org.Hs.eg.db', 'ENTREZID') %>% enrichplot::pairwise_termsim(.) %>% treeplot(.) + scale_color_viridis_c() + scale_fill_brewer(palette="Set1")

grid.arrange(p3, p4, ncol=1)

p3


treeplot(go_upx2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

go_upx <- setReadable(go_up, 'org.Hs.eg.db', 'ENTREZID')
geneList2 <- data.frame(subset(shr_ht_r2ad, shr_ht_r2ad$gene%in% ht_unique_up))
geneList2 <- geneList2[ , c(6, 2)]
cnetplot(go_upx, foldChange=geneList2)

mito <- c( "TOMM34", "TIMM13", "FBXO7", "SAMM50", "GZMB",
           "GDAP1", "PDCD5", "TIMM50", "GRPEL1", "MTCH2", 
           "RHOU", "ROMO1", "TOMM40", "TIMM10B", "AKT1",
           "TIMM8B", "TOMM5", "LEPROT")

pheatmap(assay(vsd_r2add)[mito, ], scale="row")
View(data.frame(shr_ht_r2ad_sig))

### Checking expression of some genes:

dds_r2adb <- dds_r2ad
rownames(dds_r2adb) <- mcols(dds_r2ad)$gene_name

countsplot <- function(x){
  d <- plotCounts(dds_r2adb, gene=x, intgroup="genotype", returnData = T, normalized=T)
  ggplot(d, aes(x = genotype, y = count, fill = genotype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(size=3, alpha=0.5, position=position_jitterdodge(jitter.width=0.3, dodge.width = 2))+
    theme_classic() +
    ggtitle(x) +
    xlab("Genotype")+
    ylab("Normalized counts")+
    labs(fill = "Genotype")+
    scale_fill_brewer(palette="Set1")+
    theme(plot.title = element_text(hjust = 0.5))}

countsplot("TTK")
CCgenes <- list("TTK", "BUB1B", "CDC14A", "BUB1", "PCNA", "CDC25B", "RB1", "CDC7", "MCM6", "FZR1", "ANAPC5", "ABL1", "CDK4", "TGFB1")
mat_r2ad <- assay(vsd_r2adb)
mat_r2adb <- mat_r2ad
rownames(mat_r2adb) <- mcols(dds_r2ad)$gene_name

vstcountplot <- function(x){mat_r2adb %>%
    subset(. , rownames(.)==x) %>%
    as.data.frame(.)%>%
    gather()%>%
    mutate(Genotype = factor(c(rep("WT", 6), rep("Het", 5)), levels=c("WT", "Het"), ordered=T)) %>%
    ggplot(aes(x=Genotype, y=value))+
    geom_boxplot(aes(fill=Genotype), alpha=0.8,outlier.shape = NA)+
    geom_point(color="black", size=3, alpha=0.5)+
    theme_classic()+
    scale_fill_brewer(palette="Set1")+
    labs(title=x, x="", y="Normalized counts")}

genesup <- ht_unique_res_hbm %>%
  subset(. , log2FoldChange>0.58) %>%
  subset(. , lfcSE<0.2)%>% as.data.frame(.) %>% .$gene

genesdown <- ht_unique_res_hbm %>%
  subset(. , log2FoldChange<c(-0.58)) %>%
  subset(. , lfcSE<0.2)%>% as.data.frame(.) %>% .$gene

View(data.frame(genesup)) ## See if there are any duplicates as these mess up the vstcountsplot function

genesupplots <- lapply(genesup[-c(19, 24, 36)], function(x){
  vstcountplot(x)})
grid.arrange(grobs=genesupplots, ncol=8)

genesdownplots <- lapply(genesdown, function(x){
  vstcountplot(x)})
grid.arrange(grobs=genesdownplots, ncol=3)

vstcountplot("PAX6")

gg <- c("PAX6", "CCND1", "NEUROD1", "CNTN2")

ggs <- lapply(gg, function(x){
  vstcountplot(x)})
grid.arrange(grobs=ggs, ncol=2)

# GSVA
# View(msigdbr_collections())

BPs <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
BP_list <- split(BPs$entrez_gene,
                 BPs$gs_name)
CCs <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "CC")
CC_list <- split(CCs$entrez_gene,
                 CCs$gs_name)
MFs <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")
MF_list <- split(MFs$entrez_gene,
                 MFs$gs_name)
Hs <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
H_list <- split(Hs$entrez_gene,
                Hs$gs_name)
BIOCs <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA")
BIOC_list <- split(BIOCs$entrez_gene,
                   BIOCs$gs_name)
KEGGs <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
KEGG_list <- split(KEGGs$entrez_gene,
                   KEGGs$gs_name)
REACTs <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
REACT_list <- split(REACTs$entrez_gene,
                    REACTs$gs_name)
WIKIs <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS")
WIKI_list <- split(WIKIs$entrez_gene,
                   WIKIs$gs_name)
PIDs <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:PID")
PID_list <- split(PIDs$entrez_gene,
                  PIDs$gs_name)
CGPs <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
CGP_list <- split(CGPs$entrez_gene,
                  CGPs$gs_name)
C4s <- msigdbr::msigdbr(species = "Homo sapiens", category = "C4")
C4_list <- split(C4s$entrez_gene,
                 C4s$gs_name)
C6s <- msigdbr::msigdbr(species = "Homo sapiens", category = "C6")
C6_list <- split(C6s$entrez_gene,
                 C6s$gs_name)

# Retrieve the normalized data from the `DESeqDataSet`
ht0_df <- assay(vsd_r2adb) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ensembl_id")

ht0_mapped <- data.frame(
  "entrez_id" = mapIds(
    org.Hs.eg.db,
    keys = ht0_df$ensembl_id,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first")) %>%  
  dplyr::filter(!is.na(entrez_id)) %>%
  tibble::rownames_to_column("Ensembl") %>%
  dplyr::inner_join(ht0_df, by = c("Ensembl" = "ensembl_id"))

gene_means <- rowMeans(ht0_mapped %>% dplyr::select(-Ensembl, -entrez_id))

# Let's add this as a column in our `ht0_mapped`.
ht0_mapped <- ht0_mapped %>%
  dplyr::mutate(gene_means) %>%
  dplyr::select(Ensembl, entrez_id, gene_means, dplyr::everything())

ht0_filt <- ht0_mapped %>%
  dplyr::arrange(dplyr::desc(gene_means)) %>%
  dplyr::distinct(entrez_id, .keep_all = TRUE)

sum(duplicated(ht0_filt$entrez_id))

ht0_filt_mat <- ht0_filt %>%
  dplyr::select(-Ensembl, -gene_means) %>%
  tibble::column_to_rownames("entrez_id") %>%  as.matrix()

ht0_gsva_BP <- gsva(ht0_filt_mat, BP_list, method = "gsva", kcdf = "Gaussian",min.sz = 15, max.sz = 500, mx.diff = TRUE)
ht0_gsva_CC <- gsva(ht0_filt_mat, CC_list, method = "gsva", kcdf = "Gaussian",min.sz = 15, max.sz = 500, mx.diff = TRUE)
ht0_gsva_MF <- gsva(ht0_filt_mat, MF_list, method = "gsva", kcdf = "Gaussian",min.sz = 15, max.sz = 500, mx.diff = TRUE)
ht0_gsva_H <- gsva(ht0_filt_mat, H_list, method = "gsva", kcdf = "Gaussian",min.sz = 15, max.sz = 500, mx.diff = TRUE)
ht0_gsva_BIOC <- gsva(ht0_filt_mat, BIOC_list, method = "gsva", kcdf = "Gaussian",min.sz = 15, max.sz = 500, mx.diff = TRUE)
ht0_gsva_KEGG <- gsva(ht0_filt_mat, KEGG_list, method = "gsva", kcdf = "Gaussian",min.sz = 15, max.sz = 500, mx.diff = TRUE)
ht0_gsva_REACT <- gsva(ht0_filt_mat, REACT_list, method = "gsva", kcdf = "Gaussian",min.sz = 15, max.sz = 500, mx.diff = TRUE)
ht0_gsva_WIKI <- gsva(ht0_filt_mat, WIKI_list, method = "gsva", kcdf = "Gaussian",min.sz = 15, max.sz = 500, mx.diff = TRUE)
ht0_gsva_PID <- gsva(ht0_filt_mat, PID_list, method = "gsva", kcdf = "Gaussian",min.sz = 15, max.sz = 500, mx.diff = TRUE)
ht0_gsva_CGP <- gsva(ht0_filt_mat, CGP_list, method = "gsva", kcdf = "Gaussian",min.sz = 15, max.sz = 500, mx.diff = TRUE)
ht0_gsva_C4 <- gsva(ht0_filt_mat, C4_list, method = "gsva", kcdf = "Gaussian",min.sz = 15, max.sz = 500, mx.diff = TRUE)
ht0_gsva_C6 <- gsva(ht0_filt_mat, C6_list, method = "gsva", kcdf = "Gaussian",min.sz = 15, max.sz = 500, mx.diff = TRUE)

#Setting up mm
mod <- model.matrix(~ factor(dds_r2ad$genotype))
colnames(mod) <- c("WT", "Het")

# DE pathways

#BP
fit_ht0_BP <- lmFit(ht0_gsva_BP, mod) %>%
  eBayes()
fit_ht0_CC <- lmFit(ht0_gsva_CC, mod) %>%
  eBayes()
fit_ht0_MF <- lmFit(ht0_gsva_MF, mod) %>%
  eBayes()
fit_ht0_H <- lmFit(ht0_gsva_H, mod) %>%
  eBayes()
fit_ht0_BIOC <- lmFit(ht0_gsva_BIOC, mod) %>%
  eBayes()
fit_ht0_KEGG <- lmFit(ht0_gsva_KEGG, mod) %>%
  eBayes()
fit_ht0_REACT <- lmFit(ht0_gsva_REACT, mod) %>%
  eBayes()
fit_ht0_WIKI <- lmFit(ht0_gsva_WIKI, mod) %>%
  eBayes()
fit_ht0_PID <- lmFit(ht0_gsva_PID, mod) %>%
  eBayes()
fit_ht0_CGP <- lmFit(ht0_gsva_CGP, mod) %>%
  eBayes()
fit_ht0_C4 <- lmFit(ht0_gsva_C4, mod) %>%
  eBayes()
fit_ht0_C6 <- lmFit(ht0_gsva_C6, mod) %>%
  eBayes()

# How many pathways p<0.05
decideTests(fit_ht0_BP, p.value=0.05) %>%
  summary() # het 10 down, 30 up / hom 286 down, 198 up
decideTests(fit_ht0_CC, p.value=0.05) %>%
  summary() # het 8 down 12 up / hom 31 down, 13 up
decideTests(fit_ht0_MF, p.value=0.05) %>%
  summary() # het 4 down 4 up / hom 45 down 35 up
decideTests(fit_ht0_H, p.value=0.05) %>%
  summary() # het 1 down 3 up / hom 1 down 1 up
decideTests(fit_ht0_BIOC, p.value=0.05) %>%
  summary() # het 0 down 1 up / hom 24 down 2 up
decideTests(fit_ht0_KEGG, p.value=0.05) %>%
  summary() # het 0 down 0 up / hom 10 down 15 up
decideTests(fit_ht0_REACT, p.value=0.05) %>%
  summary() # het 3 down 9 up / hom 36 down 43 up
decideTests(fit_ht0_WIKI, p.value=0.05) %>%
  summary() # het 3 down 7 up / hom 32 down 24 up
decideTests(fit_ht0_PID, p.value=0.05) %>%
  summary() # het 0 down 4 up / hom 30 down 10 up
decideTests(fit_ht0_CGP, p.value=0.05) %>%
  summary() # het 13 down 101 up / hom 91 down 150 up
decideTests(fit_ht0_C4, p.value=0.05) %>%
  summary() # het 1 down 27 up / hom 20 down 46 up
decideTests(fit_ht0_C6, p.value=0.05) %>%
  summary() # het 0 down 0 up / hom 5 down 3 up

#BP
BP_ht0_sig <- topTable(fit_ht0_BP, coef=2, n=Inf) %>%
  subset(adj.P.Val <= 0.05)
BP_ht0_t10 <- rbind(na.omit(BP_ht0_sig[order(BP_ht0_sig$logFC)[1:20],]),
                   na.omit(BP_ht0_sig[order(-BP_ht0_sig$logFC)[1:20],]))
#CC
CC_ht0_sig <- topTable(fit_ht0_CC, coef=2, n=Inf) %>%
  subset(adj.P.Val <= 0.05)
CC_ht0_t10 <- rbind(na.omit(CC_ht0_sig[order(CC_ht0_sig$logFC)[1:20],]),
                   na.omit(CC_ht0_sig[order(-CC_ht0_sig$logFC)[1:20],]))

#MF
MF_ht0_sig <- topTable(fit_ht0_MF, coef=2, n=Inf) %>%
  subset(adj.P.Val <= 0.05)
MF_ht0_t10 <- rbind(na.omit(MF_ht0_sig[order(MF_ht0_sig$logFC)[1:20],]),
                   na.omit(MF_ht0_sig[order(-MF_ht0_sig$logFC)[1:20],]))
#H
H_ht0_sig <- topTable(fit_ht0_H, coef=2, n=Inf) %>%
  subset(adj.P.Val <= 0.05)
H_ht0_t10 <- rbind(na.omit(H_ht0_sig[order(H_ht0_sig$logFC)[1:3],]),
                  na.omit(H_ht0_sig[order(-H_ht0_sig$logFC)[1:4],]))
#BIOC
BIOC_ht0_sig <- topTable(fit_ht0_BIOC, coef=2, n=Inf) %>%
  subset(adj.P.Val <= 0.05)
BIOC_ht0_t10 <- rbind(na.omit(BIOC_ht0_sig[order(BIOC_ht0_sig$logFC)[1:9],]),
                     na.omit(BIOC_ht0_sig[order(-BIOC_ht0_sig$logFC)[1:2],]))

#KEGG
KEGG_ht0_sig <- topTable(fit_ht0_KEGG, coef=2, n=Inf) %>%
  subset(adj.P.Val <= 0.05)
KEGG_ht0_t10 <- rbind(na.omit(KEGG_ht0_sig[order(KEGG_ht0_sig$logFC)[1:10],]),
                     na.omit(KEGG_ht0_sig[order(-KEGG_ht0_sig$logFC)[1:9],]))
#REACT
REACT_ht0_sig <- topTable(fit_ht0_REACT, coef=2, n=Inf) %>%
  subset(adj.P.Val <= 0.05)
REACT_ht0_t10 <- rbind(na.omit(REACT_ht0_sig[order(REACT_ht0_sig$logFC)[1:10],]),
                      na.omit(REACT_ht0_sig[order(-REACT_ht0_sig$logFC)[1:10],]))
#WIKI
WIKI_ht0_sig <- topTable(fit_ht0_WIKI, coef=2, n=Inf) %>%
  subset(adj.P.Val <= 0.05)
WIKI_ht0_t10 <- rbind(na.omit(WIKI_ht0_sig[order(WIKI_ht0_sig$logFC)[1:20],]),
                     na.omit(WIKI_ht0_sig[order(-WIKI_ht0_sig$logFC)[1:20],]))

#PID
PID_ht0_sig <- topTable(fit_ht0_PID, coef=2, n=Inf) %>%
  subset(adj.P.Val <= 0.05)
PID_ht0_t10 <- rbind(na.omit(PID_ht0_sig[order(PID_ht0_sig$logFC)[1:9],]),
                    na.omit(PID_ht0_sig[order(-PID_ht0_sig$logFC)[1:6],]))

#CGP
CGP_ht0_sig <- topTable(fit_ht0_CGP, coef=2, n=Inf) %>%
  subset(adj.P.Val <= 0.05)
CGP_ht0_t10 <- rbind(na.omit(CGP_ht0_sig[order(CGP_ht0_sig$logFC)[1:20],]),
                    na.omit(CGP_ht0_sig[order(-CGP_ht0_sig$logFC)[1:20],]))

#C4
C4_ht0_sig <- topTable(fit_ht0_C4, coef=2, n=Inf) %>%
  subset(adj.P.Val <= 0.05)
C4_ht0_t10 <- rbind(na.omit(C4_ht0_sig[order(C4_ht0_sig$logFC)[1:20],]),
                   na.omit(C4_ht0_sig[order(-C4_ht0_sig$logFC)[1:20],]))

#C6
C6_ht0_sig <- topTable(fit_ht0_C6, coef=1, n=Inf) %>%
  subset(adj.P.Val <= 0.05)
C6_ht0_t10 <- rbind(na.omit(C6_ht0_sig[order(C6_ht0_sig$logFC)[1:88],]))

DEpw_ht <- rbind(BP_ht0_sig, CC_ht0_sig, MF_ht0_sig,
                  H_ht0_sig, BIOC_ht0_sig, KEGG_ht0_sig, REACT_ht0_sig,
                  WIKI_ht0_sig, PID_ht0_sig, CGP_ht0_sig, C4_ht0_sig, C6_ht0_sig)
DEpw_ht$collection <- c(rep("GO:BP", nrow(BP_ht0_sig)), rep("GO:CC", nrow(CC_ht0_sig)), rep("GO:MF", nrow(MF_ht0_sig)),
                        rep("HALLMARK", nrow(H_ht0_sig)), rep("BIOCARTA", nrow(BIOC_ht0_sig)), rep("KEGG", nrow(KEGG_ht0_sig)), 
                        rep("REACTOME", nrow(REACT_ht0_sig)), rep("WIKIPATHWAYS", nrow(WIKI_ht0_sig)), rep("PID", nrow(PID_ht0_sig)), 
                        rep("CGP", nrow(CGP_ht0_sig)), rep("C4", nrow(C4_ht0_sig)), rep("C6", nrow(C6_ht0_sig)))

all_gene_sets <- rbind(BPs, CCs, MFs, Hs, BIOCs, KEGGs, REACTs, WIKIs, PIDs, CGPs, C4s, C6s)

write.csv(DEpw_ht, file="DEpws")
write.csv(all_gene_sets, file="all_gene_sets")

View(DEpw_ht)

DEpw_ht0 <- rbind(BP_ht0_t10, CC_ht0_t10, MF_ht0_t10,
                 H_ht0_t10, BIOC_ht0_t10, KEGG_ht0_t10, REACT_ht0_t10,
                 WIKI_ht0_t10, PID_ht0_t10, CGP_ht0_t10, C4_ht0_t10, C6_ht0_t10)
DEpw_ht0$Geneset <- c(rep("BP", nrow(BP_ht0_t10)), rep("CC", nrow(CC_ht0_t10)), rep("MF", nrow(MF_ht0_t10)),
                     rep("H", nrow(H_ht0_t10)), rep("BIOC", nrow(BIOC_ht0_t10)), rep("KEGG", nrow(KEGG_ht0_t10)),
                     rep("REACT", nrow(REACT_ht0_t10)), rep("WIKI", nrow(WIKI_ht0_t10)),
                     rep("PID", nrow(PID_ht0_t10)), rep("CGP", nrow(CGP_ht0_t10)), 
                     rep("C4", nrow(C4_ht0_t10)), rep("C6", nrow(C6_ht0_t10)))


BPplot_ht0 <- ggplot(subset(DEpw_ht0, Geneset=="BP"), aes(y=reorder(rownames(subset(DEpw_ht0, Geneset=="BP")), logFC), x=logFC)) + 
  geom_col(aes(fill=adj.P.Val)) +
  ylab("Biological pathway")
CCplot_ht0 <- ggplot(subset(DEpw_ht0, Geneset=="CC"), aes(y=reorder(rownames(subset(DEpw_ht0, Geneset=="CC")), logFC), x=logFC)) + 
  geom_col(aes(fill=adj.P.Val)) +
  ylab("Cellular Component")
MFplot_ht0 <- ggplot(subset(DEpw_ht0, Geneset=="MF"), aes(y=reorder(rownames(subset(DEpw_ht0, Geneset=="MF")), logFC), x=logFC)) + 
  geom_col(aes(fill=adj.P.Val)) +
  ylab("Molecular Function")
Hplot_ht0 <- ggplot(subset(DEpw_ht0, Geneset=="H"), aes(y=reorder(rownames(subset(DEpw_ht0, Geneset=="H")), logFC), x=logFC)) + 
  geom_col(aes(fill=adj.P.Val)) +
  ylab("Hallmark")
BIOCplot_ht0 <- ggplot(subset(DEpw_ht0, Geneset=="BIOC"), aes(y=reorder(rownames(subset(DEpw_ht0, Geneset=="BIOC")), logFC), x=logFC)) + 
  geom_col(aes(fill=adj.P.Val)) +
  ylab("BIOCARTA")
KEGGplot_ht0 <- ggplot(subset(DEpw_ht0, Geneset=="KEGG"), aes(y=reorder(rownames(subset(DEpw_ht0, Geneset=="KEGG")), logFC), x=logFC)) + 
  geom_col(aes(fill=adj.P.Val)) +
  ylab("KEGG")
REACTplot_ht0 <- ggplot(subset(DEpw_ht0, Geneset=="REACT"), aes(y=reorder(rownames(subset(DEpw_ht0, Geneset=="REACT")), logFC), x=logFC)) + 
  geom_col(aes(fill=adj.P.Val)) +
  ylab("REACTOME")
WIKIplot_ht0 <- ggplot(subset(DEpw_ht0, Geneset=="WIKI"), aes(y=reorder(rownames(subset(DEpw_ht0, Geneset=="WIKI")), logFC), x=logFC)) + 
  geom_col(aes(fill=adj.P.Val)) +
  ylab("WIKIathways")
PIDplot_ht0 <- ggplot(subset(DEpw_ht0, Geneset=="PID"), aes(y=reorder(rownames(subset(DEpw_ht0, Geneset=="PID")), logFC), x=logFC)) + 
  geom_col(aes(fill=adj.P.Val)) +
  ylab("PID")
CGPplot_ht0 <- ggplot(subset(DEpw_ht0, Geneset=="CGP"), aes(y=reorder(rownames(subset(DEpw_ht0, Geneset=="CGP")), logFC), x=logFC)) + 
  geom_col(aes(fill=adj.P.Val)) +
  ylab("Cancer computational gene sets")
C4plot_ht0 <- ggplot(subset(DEpw_ht0, Geneset=="C4"), aes(y=reorder(rownames(subset(DEpw_ht0, Geneset=="C4")), logFC), x=logFC)) + 
  geom_col(aes(fill=adj.P.Val)) +
  ylab("C4")
C6plot_ht0 <- ggplot(subset(DEpw_ht0, Geneset=="C6"), aes(y=reorder(rownames(subset(DEpw_ht0, Geneset=="C6")), logFC), x=logFC)) + 
  geom_col(aes(fill=adj.P.Val)) +
  ylab("C6")
View(BP_ht0_t10)
BPplot_ht0
CCplot_ht0
MFplot_ht0
Hplot_ht0
BIOCplot_ht0
KEGGplot_ht0
REACTplot_ht0
WIKIplot_ht0
PIDplot_ht0
CGPplot_ht0
C4plot_ht0
C6plot_ht0

save.image("/Volumes/Partitie4/RNAseq_analysis_2/RNAseqR2_2.RData")

rel_PW <- c("GOBP_NUCLEOTIDE_EXCISION_REPAIR_DNA_GAP_FILLING", "GOBP_KINETOCHORE_ORGANIZATION", 
            "GOBP_POSITIVE_REGULATION_OF_NEURAL_PRECURSOR_CELL_PROLIFERATION", "GOBP_CELL_PROLIFERATION_IN_HINDBRAIN", 
            "GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE_G2_M_PHASE_TRANSITION", "GOBP_PROTEIN_IMPORT_INTO_MITOCHONDRIAL_MATRIX", 
            "GOBP_PSEUDOURIDINE_SYNTHESIS", "GOBP_CELL_JUNCTION_DISASSEMBLY", 
            "GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS", 
            "GOBP_NEGATIVE_REGULATION_OF_CELL_MORPHOGENESIS_INVOLVED_IN_DIFFERENTIATION")
rel_sub <- subset(BP_ht0_sig, rownames(BP_ht0_sig)%in%rel_PW)
rownames(rel_sub) <- case_when(rownames(rel_sub) =="GOBP_KINETOCHORE_ORGANIZATION" ~"Kinetochore organisation", 
                               rownames(rel_sub) == "GOBP_NUCLEOTIDE_EXCISION_REPAIR_DNA_GAP_FILLING" ~ "Nucleotide excisino repair DNA gap filling", 
                               rownames(rel_sub) =="GOBP_CELL_PROLIFERATION_IN_HINDBRAIN" ~ "Cell proliferation in hindbrain", 
                               rownames(rel_sub) =="GOBP_POSITIVE_REGULATION_OF_CELL_CYCLE_G2_M_PHASE_TRANSITION" ~ "Positive regulation of cell cycle G2-M phase transition",
                               rownames(rel_sub) =="GOBP_POSITIVE_REGULATION_OF_NEURAL_PRECURSOR_CELL_PROLIFERATION" ~ "Positive regulation of neural precursor cell proliferation",
                               rownames(rel_sub) =="GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS" ~ "Negative regulation of Extrinsic apopotitic signalling pathway via death domain receptors",
                               rownames(rel_sub) =="GOBP_CELL_JUNCTION_DISASSEMBLY" ~ "Cell junction disassembly",
                               rownames(rel_sub) =="GOBP_PSEUDOURIDINE_SYNTHESIS" ~ "Pseudouridine synthesis", 
                               rownames(rel_sub) =="GOBP_PROTEIN_IMPORT_INTO_MITOCHONDRIAL_MATRIX" ~ "Protein import into mitochondrial matrix", 
                               rownames(rel_sub) =="GOBP_NEGATIVE_REGULATION_OF_CELL_MORPHOGENESIS_INVOLVED_IN_DIFFERENTIATION" ~ "Negative regulation of cell morphogenesis involved in differentiation")
rel_ht0 <- ggplot(rel_sub, aes(y=reorder(rownames(rel_sub), logFC), x=logFC)) + 
  geom_col(aes(fill=adj.P.Val)) +
  ylab("")+ theme(axis.text = element_text(size = 10))
rel_ht0
#### Pathways of interest:
GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS

## Exploring the genes:
shr_ht_r2ad$gene <- mcols(dds_r2ad)$gene_name
View(data.frame(shr_ht_r2ad))

PWtotl <- rbind(BPs, CCs, MFs, Hs, BIOCs, KEGGs, REACTs, WIKIs, PIDs, C4s, C6s, CGPs)

KCC <- PWtotl[PWtotl$gs_name %in% "KEGG_CELL_CYCLE", ]
KCCsig <- subset(shr_ht_r2ad_sig, rownames(shr_ht_r2ad_sig)%in%KCC$ensembl_gene)
CPH <- PWtotl[PWtotl$gs_name %in% "GOBP_CELL_PROLIFERATION_IN_HINDBRAIN", ]
CPHsig <- subset(shr_ht_r2ad_sig, rownames(shr_ht_r2ad_sig)%in%CPH$ensembl_gene)
NAP <- PWtotl[PWtotl$gs_name %in% "GOBP_NEGATIVE_REGULATION_OF_EXTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_VIA_DEATH_DOMAIN_RECEPTORS", ]
NAPsig <- subset(shr_ht_r2ad_sig, rownames(shr_ht_r2ad_sig)%in%NAP$ensembl_gene)

vsd_r2add <- vsd_r2adb
rownames(vsd_r2add) <- mcols(dds_r2ad)$gene_name
mat_r2ad <- subset(assay(vsd_r2add), rowVars(assay(vsd_r2add))!=0)
cbdev <- c("BMP6", "BMP7", "ATOH1", "MEIS1", "PAX6", "ZIC1", "ZIC3", "SHH", "NOTCH2", "JAG1", "ATOH1", "MYCN", "BMP4", "WNT3", "APC2", "APC11", "NEUROD1", "ZIC2", "SEMA6A", "CNTN2", "PARD3", "JAM3", "ASTN1", "PARD6")
cycling <- c("ZIC1", "PTCH1", "GLI1", "MKI67", "NEUROD1", "CNTN2", "SEMA6A", "RBFOX3")

pheatmap(mat_r2ad[rownames(mat_r2ad)%in%cbdev ], scale="row",
         cluster_rows = T, cluster_cols = F,
         labels_col=cano, annotation_col=ca)


pheatmap(mat_r2ad[rownames(mat_r2ad)%in%cbdev, ], scale="row", cluster_rows = T, cluster_cols = F, labels_col=cano, annotation_col=ca)


PWsig <- subset(shr_ht_r2ad_sig, rownames(shr_ht_r2ad_sig)%in%PW$ensembl_gene)
PWall <- subset(shr_ht_r2ad, rownames(shr_ht_r2ad)%in%PW$ensembl_gene)
pheatmap(mat_r2ad[rownames(mat_r2ad)%in%rownames(PWsig), ], scale="row", 
         labels_row=PWsig$gene, labels_col=ca2, annotation_col=ca,
         main="KEGG", show_rownames=F)

PWht_unique <- subset(PWsig, rownames(PWsig) %in% ht_unique)
View(data.frame(PWht_unique))

### Volcanoplot Cellcycle genes
CCgenes <- list("MKI67", "BUB1", "PCNA",  "RB1", "CDC7", "MCM6")


ggplot(data.frame(shr_ht_r2ad), aes(y=baseMean, x=log2FoldChange, color=gene_change))+
  geom_point()+
  scale_y_log10(breaks=c(100, 250, 500, 1000, 1500, 2500, 5000, 7500, 10000), limits=c(250, 35000))+
  scale_x_continuous(limits=c(-4, 6))+
  geom_label_repel(data = labels5, aes(label = gene, color=gene_change), max.overlaps = Inf)

### Volcano plots etc:
test2 <- data.frame(shr_ht_r2ad)  %>%
  mutate(gene_change = case_when(log2FoldChange >= 1 & baseMean > 500 ~ "Up high",
                                 log2FoldChange <= -1 & baseMean > 500 ~ "Down low",
                                 TRUE ~ "ns"))  %>%
  mutate(gene = mcols(dds_r2ad)$gene_name)

labels3 <- subset(test2, test2$gene_change!="ns")
ggplot(test2, aes(y=baseMean, x=log2FoldChange, color=gene_change))+
  geom_point()+
  scale_y_log10(breaks=c(100, 250, 500, 1000, 1500, 2500, 5000, 7500, 10000), limits=c(250, 35000))+
  scale_x_continuous(limits=c(-4, 6))+
  geom_label_repel(data = labels3, aes(label = gene, color=gene_change), max.overlaps = Inf)

test3 <- data.frame(shr_ht_r2ad)  %>%
  mutate(gene_change = case_when(log2FoldChange >= 0.585 & baseMean > 300 ~ "Up high",
                                 log2FoldChange <= -0.585 & baseMean > 300 ~ "Down low",
                                 TRUE ~ "ns"))  %>%
  mutate(gene = mcols(dds_r2ad)$gene_name)

labels4 <- subset(test3, rownames(test3)%in%rownames(ht_unique_res_a_hbm))
labels5 <- subset(labels4, labels4$gene_change!="ns")
ggplot(test3, aes(y=baseMean, x=log2FoldChange, color=gene_change))+
  geom_point()+
  scale_y_log10(breaks=c(100, 250, 500, 1000, 1500, 2500, 5000, 7500, 10000), limits=c(250, 35000))+
  scale_x_continuous(limits=c(-4, 6))+
  geom_label_repel(data = labels5, aes(label = gene, color=gene_change), max.overlaps = Inf)

### Comparison with recently published snRNAseq from ptch1het mice:
ptch1up <- read_xlsx("/Users/max/Documents/ptch1up.xlsx")
human_symbols = convert_mouse_to_human_symbols(ptch1up$gene)
df <- data.frame(ptch1up$gene)
df$human <- human_symbols
colnames(df) <- c("mouse", "human")
ptch1up_b <- inner_join(ptch1up, df, by=c("gene"="mouse"))
ptch1up_c <- ptch1up_b[,c(8, 3,7)]
ptch1up_d <- subset(ptch1up_c, ptch1up_c$cell_type=="Dividing GNPs-like tumor cells")

ptch1up_genes <- shr_ht_r2ad_sig %>%
  subset(., gene%in%ptch1up_d$human) %>% .$gene

ints <- shr_ht_r2ad_sig %>%
  subset(., gene%in%ptch1up_d$human)

vsd_r2add <- vsd_r2adb
rownames(vsd_r2add) <- mcols(dds_r2ad)$gene_name
pheatmap(assay(vsd_r2add)[ints$gene, ], scale="row", 
         show_rownames = F, labels_col=ca2, annotation_col=ca,
         main="Genes upregulated in PTCH1+/- mouse 
         dividing GNPs-like tumor cells")

ints_unique <- ht_unique_res %>%
  subset(., gene%in%ptch1up_d$human) %>%
  subset(., log2FoldChange>0)

pheatmap(assay(vsd_r2add)[ints_unique$gene, ], scale="row", 
         show_rownames = T, labels_col=ca2, annotation_col=ca,
         main="Genes upregulated in PTCH1+/- mouse 
         dividing GNPs-like tumor cells")

pheatmap(assay(vsd_r2add)[rownames(vsd_r2add)%in%labels, ], scale="row")

# #####
# labels <- c("PTCH1", "GLI1", "GLI2", "NEUROD1", "CNTN2")
labels <- c("SHH", "PTCH1", "SMO", "GLI1", "HIF1A", "ZEB1", "ITGB1", "PARD3", "CXCL3", "NEUROD1", "ZIC1", "ZIC2", "PAX6") 
labels2 <- c( "CNTN2", "NEUROD2", "ZIC1", "ZIC2", "SEMA6A", "PLXNA2", "CNTN1", "PARD3", "PAX6")
labels3 <- c("ASTN1", "ASTN2", "CDH2", "PAR6", "FZR1")
labels4 <- c("PAX6", "PLXNA2", "RBFOX3")
# labels <- c("SSTR1", "SSTr2ad", "SSTR3", "SSTR4", "SSTR5", "ADCYAP1R1", "VIPR1", "CNTN2", "GABRA1", "GABRA2", "GABRA3", "GABRA6", "GABRB2", "GABRB3", "GABRD", "GABRG1", "GABRG2", "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C" , "RELN")
# labels <- c("CENPE", "TOP2A", "HMMR", "DLGAP5", "HDAC9")
# labels <- c("SPP1", "NID2", "UNC5C", "PTN", "COL4A1", "PAX6", "CD47", "PTPRZ1", "GJA1", "PLA2G7")
# labels <- c("BAG3", "ATF3", "HSPH1", "DNAJB1", "DNAJB2", "MTX2")
# labels <- c("FOSB", "FOXF2", "KLF4", "CSRP1", "ZBTB20", "EGr2ad", "SOX2", "NEUROD1", "NDN", "ZIC3")
# labels <- c("PAX6", "MEIS1", "ZIC1", "ZIC2", "BARHL1", "ATOH1", "NEUROD1", "UNC5C", "CXCR4", "SEMA6A", "DCX")
# labels  <- c("ATOH1", "MYCN", "CCND2", "GLI1", "NOTCH2", "MKI67", "PHH3")
# labels <- c("NEUROD1", "DCX", "CNTN2", "NSCL1", "ZIC2", "FZR1", "CDKN1B")
 labels <- c("SEMA6A", "PLXNA2", "EN2", "ASTN1", "ASTN2", "PARD6A", "CDK5")
 labels <- c("NEUROD1", "CNTN2", "SPP1", "NID2", "UNC5C", "PAX6", "CD47", "PTPRZ1", "PCNA", "PLK1", "TOP2A", "CCND1", "CCNB1", "MKI67",
             "CENPE", "HDAC9")
 

# 
# 
test <- data.frame(shr_ht_r2ad)  %>%
  mutate(gene_type= case_when(log2FoldChange>=0.5849625 & padj <=0.05 ~ "Up", 
                              log2FoldChange<= c(-0.5849625) & padj<=0.05 ~ "Down", 
                              TRUE ~ "ns")) %>%
  mutate(gene = mcols(dds_r2ad)$gene_name)

SHH <- c("SHH", "PTCH1", "PTCH2", "GLI1", "GLI2", "GLI3", "SMO")
SHHlab <- subset(test, test$gene %in% SHH)
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 

SHH_vol <- ggplot(data=test, aes(x = log2FoldChange,
                      y = -log10(padj)))+
  geom_point(alpha=0.3, aes(color=gene_type)) +
  geom_point(data = SHHlab,
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
  scale_x_continuous(limits=c(-5,5)) +
  geom_label_repel(data = SHHlab,   
                   aes(label = gene),
                   max.overlaps = Inf) +
  scale_color_manual(values = c("#26b3ff", "grey", "#ffad73"))+
  labs(title = "SHH pathway genes in PTCH1+/- vs WT organoids",
       x = "log2(fold change)",
       y = "-log10(BH adjusted P-value)",
       colour = "Expression \nchange") +
  theme_classic()

SHH_H <- pheatmap(mat_r2adb[SHH, ], scale="row", annotation_col=ca, cellwidth=12, cellheight=12, treeheight_col=4, treeheight_row=4, border_color="white", show_colnames = F, legend=F)


SHH <- grid.arrange(SHH_vol, SHH_H[[4]], ncol=2)


topRLsvz <- data.frame(RLSVZ_sig[order(RLSVZ_sig$log2FoldChange), ])[1:100, ]$gene
topRLvz <- data.frame(RLVZ_sig[order(RLVZ_sig$log2FoldChange), ])[1:100, ]$gene
RLsvz <- intersect(RLSVZ_up$gene,  shr_ht_r2ad_up$gene)
RLvz <- intersect(RLVZ_up$gene,  shr_ht_r2ad_up$gene)
Egl <- intersect(EGL_up$gene, shr_ht_r2ad_up$gene)
RLsvzlab <- subset(test, test$gene %in% RLsvz)
RLvzlab <- subset(test, test$gene %in% RLvz)
EGLlab <- subset(test, test$gene %in% Egl)
RLlabs <- rbind(RLsvzlab, RLvzlab, EGLlab)

RLsvzdwn <- intersect(RLSVZ_dwn$gene,  shr_ht_r2ad_down$gene)
RLvzdwn <- intersect(RLVZ_dwn$gene,  shr_ht_r2ad_down$gene)
Egldwn <- intersect(EGL_dwn$gene, shr_ht_r2ad_down$gene)
RLsvzlabdwn <- subset(test, test$gene %in% RLsvzdwn)
RLvzlabdwn <- subset(test, test$gene %in% RLvzdwn)
EGLlabdwn <- subset(test, test$gene %in% Egldwn)
RLdwnlabs <- rbind(RLsvzlabdwn, RLvzlabdwn, EGLlabdwn)

cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
RL_mat <- rbind(RLlabs, RLdwnlabs)
RL_vol <- ggplot(data=test, aes(x = log2FoldChange,
                                 y = -log10(padj)))+
  geom_point(alpha=0.3, aes(color=gene_type)) +
  geom_point(data = RLlabs,
             size = 2,
             shape = 21,
             fill = "firebrick")  +
  geom_point(data = RLdwnlabs,
             size = 2,
             shape = 21,
             fill = "blue")  +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             alpha=0.5) +
  geom_vline(xintercept = c(-0.585, 0.585),
             linetype = "dashed",
             alpha=0.5)+
  geom_label_repel(data = RL_mat[RL_mat$gene_type=="Up"|RL_mat$gene_type=="Down", ],   
                   aes(label = gene),
                   force = 2,
                   nudge_y = 1) +
  scale_x_continuous(limits=c(-5,5)) +
  scale_color_manual(values = c("#26b3ff", "grey", "#ffad73"))+
  labs(title = "RL pathway genes in PTCH1+/- vs WT organoids",
       x = "log2(fold change)",
       y = "-log10(BH adjusted P-value)",
       colour = "Expression \nchange") +
  theme_classic()
RL_vol
tt <- rbind(RLlabs)

RL_H <- pheatmap(mat_r2adb[rownames(mat_r2adb) %in% RL_mat$gene, ], scale="row", annotation_col=ca)

RL <- grid.arrange(RL_vol, RL_H[[4]], ncol=2)


RLs <- c("PAX6",  "ATOH1", "EOMES", "TBR1", "WLS",  "SOX2","CNTN2", "NEUROD1", "NEUROD2", "ZIC1", "ZIC2", "SEMA6A", "PLXNA2", "CNTN1", "PARD3" )
RLslab <- subset(test, test$gene %in% RLs)
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 

RLs_vol <- ggplot(data=test, aes(x = log2FoldChange,
                                 y = -log10(padj)))+
  geom_point(alpha=0.3, aes(color=gene_type)) +
  geom_point(data = RLslab,
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
  scale_x_continuous(limits=c(-5,5)) +
  geom_label_repel(data = RLslab,   
                   aes(label = gene),
                   max.overlaps = Inf) +
  scale_color_manual(values = c("#26b3ff", "grey", "#ffad73"))+
  labs(title = "RL-lineage and GCP genes in PTCH1+/- vs WT organoids",
       x = "log2(fold change)",
       y = "-log10(BH adjusted P-value)",
       colour = "Expression \nchange") +
  theme_classic()

RLs_H <- pheatmap(mat_r2adb[RLs, ], scale="row", cluster_cols =F, annotation_col=ca)

EGLs <- c("PAX6", "SOX2", "WLS", "ZIC1", "ZIC2", "CNTN1", "CNTN2", "UNC5C", "NEUROD1", "NEUROD2")
ra5 <- data.frame(c(rep("Early", 5), rep("Late", 5)))
rownames(ra5) <- EGLs
colnames(ra5) <- "Stage"
pheatmap(mat_r2adb[EGLs, ], scale="row", cluster_cols =F, cluster_rows=F, annotation_col=ca, annotation_row=ra5, cellwidth=12, cellheight=12, treeheight_col=4, treeheight_row=4, border_color="white", show_colnames = F, legend=F)


RLs <- grid.arrange(RLs_vol, RLs_H[[4]], ncol=2)

CCs <- c( "PCNA", "PLK1", "TOP2A", "CCND1", "CCNB1", "MKI67", "CENPE", "CENPF", "CCNA1", "CCNA2", "CDC7", "CDCA3", "CDCA8" )
CCslab <- subset(test, test$gene %in% CCs)
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 

CCs_vol <- ggplot(data=test, aes(x = log2FoldChange,
                                 y = -log10(padj)))+
  geom_point(alpha=0.3, aes(color=gene_type)) +
  geom_point(data = CCslab,
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
  scale_x_continuous(limits=c(-5,5)) +
  geom_label_repel(data = CCslab,   
                   aes(label = gene),
                   max.overlaps = Inf) +
  scale_color_manual(values = c("#26b3ff", "grey", "#ffad73"))+
  labs(title = "Cell cycle genes in PTCH1+/- vs WT organoids",
       x = "log2(fold change)",
       y = "-log10(BH adjusted P-value)",
       colour = "Expression \nchange") +
  theme_classic()

CCs_H <- pheatmap(mat_r2adb[CCs, ], scale="row", cluster_cols =F, annotation_col=ca2, cellwidth=12, cellheight=12, treeheight_col=4, treeheight_row=4, border_color="white", show_colnames = F, legend=F)

CCs <- grid.arrange(CCs_vol, CCs_H[[4]], ncol=2)

cilios <-cilio$gene_symbol
cilioslab <- subset(test, test$gene %in% cilios)
cilioslabsig <- cilioslab[cilioslab$log2FoldChange>=0.5849625, ]
cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 

cilios_vol <- ggplot(data=test, aes(x = log2FoldChange,
                                 y = -log10(padj)))+
  geom_point(alpha=0.3, aes(color=gene_type)) +
  geom_point(data = cilioslab,
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
  scale_x_continuous(limits=c(-5,5)) +
  geom_label_repel(data = cilioslabsig,   
                   aes(label = gene),
                   max.overlaps = Inf) +
  scale_color_manual(values = c("#26b3ff", "grey", "#ffad73"))+
  labs(title = "Ciliopathy genes in PTCH1+/- vs WT organoids",
       x = "log2(fold change)",
       y = "-log10(BH adjusted P-value)",
       colour = "Expression \nchange") +
  theme_classic()

cilios_H <- pheatmap(mat_r2adb[rownames(mat_r2adb) %in% cilios, ] , scale="row", cluster_cols =F, annotation_col=ca, show_rownames = F)

cilios <- grid.arrange(cilios_vol, cilios_H[[4]], ncol=2)



# 
# dds_r2adb <- dds_r2ad
# rownames(dds_r2adb) <- mcols(dds_r2ad)$gene_name
# d <- plotCounts(dds_r2adb, gene="GLI1", intgroup="genotype", returnData = T)
# ggplot(d, aes(x = genotype, y = count, color = genotype)) + 
#   geom_point(position=position_jitter(w = 0.1,h = 0)) +
#   geom_text_repel(aes(label = rownames(d))) + 
#   theme_bw() +
#   ggtitle("MOV10") +
#   theme(plot.title = element_text(hjust = 0.5))

#### Comparison to Sam's data:
library(Seurat)
library(SeuratObject)

Sam <- readRDS("20190807vs2_v3update_newganoidscombined.rds")
Sam2 <- readRDS("cbl_integrated_cleanCC_210111.rds")
Sam
Sam <- ScaleData(Sam, vars.to.regress = "percent.mt")
DimPlot(Sam, reduction = "umap")
DimHeatmap(Sam, dims = 1, cells = 500, balanced = TRUE)

Sam_sub <- subset(Sam, features=ht_unique_res$gene)
Sam_sub_avg <- AverageExpression(Sam_sub)

counts <- data.frame(Sam_sub_avg$RNA)
counts <- subset(counts, rowVars(counts)!=0)
View(counts)
pheatmap(counts, scale="row", show_rownames=F, cluster_cols = T)
top50up <- ht_unique_res[order(ht_unique_res$log2FoldChange), ] %>% .[1:50, ]
pheatmap(counts[rownames(counts) %in%top50up$gene, ], scale="row", show_rownames=T, cluster_cols = T)

topVarGenes <- head(order(rowVars(counts), decreasing = TRUE), 50)
counts_topvar <- counts[topVarGenes,]
pheatmap(counts_topvar, scale="row", show_rownames=T, cluster_cols = T)

library(tidySingleCellExperiment)
library(scater)

Sepp <- readRDS("/Volumes/Partitie4/hum_sce.rds")
assayNames(Sepp)
normcounts(Sepp) <- log2(assay(Sepp, 'umi') + 1)
Sepp_up <- Sepp[rownames(Sepp) %in% ht_unique_up, ]
normcounts(Sepp_up) <- log2(assay(Sepp_up, 'umi') + 1)
Sepp_up
as.Seurat()

counts <- assay(Sepp_up, 'umi')
class(counts)


plotHeatmap(Sepp_up, features=rownames(Sepp_up[1:20,]), exprs_values = "normcounts")

Sepp %>%
  join_features(features=high_ht) %>%
  group_by(precisest_label) %>%
  heatmap(.feature, .cell, .abundance_counts, .scale="column")


### Heatmaps wt vs hom

FB <- c("FOXG1", "NKX2-1", "SIX3")
HB <- c("GBX2", "LMX1A", "ATOH1")
BrR <- c(FB,HB)

ra2 <- data.frame(factor(c(rep("Forebrain",  length(FB)), 
                           rep("Hindbrain", length(HB))), 
                         levels=c("Forebrain", "Hindbrain"), ordered=TRUE))
colnames(ra2)[1] <- "Region"
rownames(ra2) <- BrR
ca2 <- data.frame(vsd_r2ahm$genotype)
rownames(ca2) <- colnames(vsd_r2ahm)
colnames(ca2)[1] <- "Genotype"
mat_hom <- assay(vsd_r2ahm)
rownames(mat_hom) <- mcols(dds_r2ahm)$gene_name
pheatmap(mat_hom[BrR,], scale="row", main="Brain regions", cluster_rows=F, show_rownames=TRUE,
               cluster_cols=FALSE, annotation_row=ra2, annotation_col = ca2)

FP <- c("SHH","NTN1", "FOXA2", "NKX6-1", "NKX2-2")
RP <- c("POU4F1", "GDF7", "RSPO1", "PAX6", "WNT3A")
NeurP <- c(RP, FP)
ra3 <- data.frame(factor(c(rep("Roof plate",  length(RP)), 
                           rep("Floor plate", length(FP))), 
                         levels=c("Roof plate", "Floor plate"), ordered=TRUE))
colnames(ra3)[1] <- "Region"
rownames(ra3) <- NeurP
pheatmap(mat_hom[NeurP,], scale="row", main="Brain regions", cluster_rows=F, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_row=ra3, annotation_col = ca2)

regions <- c(NeurP, BrR)
ra4 <- rbind(ra3, ra2)
pheatmap(mat_hom[regions,], scale="row", main="Brain regions", cluster_rows=F, show_rownames=T, show_colnames = F,
         cluster_cols=FALSE, annotation_row=ra4, annotation_col = ca2, gaps_row=10,
         cellwidth = 12, cellheight = 12, border_color = "white", legend=F)



#### MB dataset exploration
MBvsd <- readRDS("/Volumes/Partitie4/RNAseq_analysis/PBTA_CBTN_MB_vsd.rds")
View(colData(MBvsd))
colannoMB <- data.frame(colData(MBvsd)$age_at_diagnosis_days, colData(MBvsd)$OS_days, colData(MBvsd)$molecular_subtype)
rownames(colannoMB) <- colnames(MBvsd)
View(colannoMB)
colnames(colannoMB) <- c("Age at diagnosis", "Overall survival", "Molecular subtype")
colannoMB$`Overall survival` <- as.numeric(colannoMB$`Overall survival`)
colannoMB$`Age at diagnosis` <- as.numeric(colannoMB$`Age at diagnosis`)

MBvsd2 <- MBvsd
rownames(MBvsd2) <- mcols(MBvsd)$gene_name
top50up <- ht_unique_res[order(-ht_unique_res$log2FoldChange), ] %>% .[1:50, ]
pheatmap(assay(MBvsd)[rownames(MBvsd)%in%rownames(top50up), ], scale="none", annotation_col = colannoMB, show_rownames = T, show_colnames = F, labels_row = top50up$gene)

top50down <- ht_unique_res[order(ht_unique_res$log2FoldChange), ] %>% .[1:50, ]
pheatmap(assay(MBvsd)[rownames(MBvsd)%in%rownames(top50down), ], scale="none", annotation_col = colannoMB, show_rownames = T, show_colnames = F, labels_row = top50down$gene)

## CBTN paediatric brain tumor exploration
CBTNdds <- readRDS("/Volumes/Partitie4/RNAseq_analysis/CBTN_u5yr_CNS/PBTA_CBTN_u5yr_dds.rds")
CBTNvsd <- readRDS("/Volumes/Partitie4/RNAseq_analysis/CBTN_u5yr_CNS/PBTA_CBTN_u5yr_vsd.rds")
View(data.frame(colData(CBTNvsd)))
colannoCBTN <- data.frame(colData(CBTNvsd)$disease_type, colData(CBTNvsd)$age_at_diagnosis_days, colData(CBTNvsd)$OS_days, colData(CBTNvsd)$molecular_subtype)
rownames(colannoCBTN) <- colnames(CBTNvsd)
View(colannoCBTN)
colnames(colannoCBTN) <- c("Disease type", "Age at diagnosis", "Overall survival", "Molecular subtype")
colannoCBTN$`Overall survival` <- as.numeric(colannoCBTN$`Overall survival`)
colannoCBTN$`Age at diagnosis` <- as.numeric(colannoCBTN$`Age at diagnosis`)
pheatmap(assay(CBTNvsd)[rownames(CBTNvsd) %in% rownames(top50down), ], 
         scale="none", annotation_col = colannoCBTN, show_rownames = T, 
         show_colnames = F, labels_row = top50down$gene)


alddsall <- readRDS("/Volumes/Partitie4/RNAseq_analysis/aldmve/aldds_all.rds")
alvsdall <- readRDS("/Volumes/Partitie4/RNAseq_analysis/aldmve/alvsd_all.rds")

al_mat <- assay(alvsdall)
rownames(al_mat) <- mcols(alddsall)$gene_name
al_mat <- data.frame(al_mat)
colnames(al_mat)
al_mat_sub1 <- al_mat[,-c(51:66)]
al_mat_sub1 <- subset(al_mat_sub1, rowVars(al_mat_sub1)!=0)
colannoALD <- data.frame(colData(alddsall)[-c(51:66),]$region, colData(alddsall)[-c(51:66),]$age_pcw)
rownames(colannoALD) <- colnames(al_mat_sub1)
colnames(colannoALD) <- c("Region", "PCW")
pheatmap(al_mat_sub1[rownames(al_mat_sub1) %in% ht_unique_res$gene, ], 
         scale="row", annotation_col = colannoALD, show_rownames = T, 
         show_colnames = F)
design(alddsall) <- ~ Age + region
alddsall <- DESeq(alddsall)
resultsNames(alddsall)
alddsall$region <- relevel(alddsall$region , ref="Bulk")
alddsall <- DESeq(alddsall)

ELGvPCL_res <- lfcShrink(alddsall, coef = "region_PCL_vs_EGL")
summary(ELGvPCL_res$padj<0.05, na.rm=TRUE) # 9786 genes
ELGvPCL_res$gene <-  mcols(alddsall)$gene_name
ELGvPCL_sig <- subset(ELGvPCL_res, ELGvPCL_res$padj<0.05)
EGLup <- subset(ELGvPCL_sig, ELGvPCL_sig$log2FoldChange<0)
View(data.frame(EGLup))
EGLuphtup <- shr_ht_r2ad_sig[rownames(shr_ht_r2ad_sig) %in%rownames(EGLup), ]$gene
View(data.frame(EGLuphtup))

pheatmap(al_mat_sub1[rownames(al_mat_sub1) %in% EGLuphtup, ], 
         scale="row", annotation_col = colannoALD, show_rownames = T, 
         show_colnames = F)

shr_ht_r2ad_up <- subset(shr_ht_r2ad_sig, shr_ht_r2ad_sig$log2FoldChange>0)
shr_ht_r2ad_down <- subset(shr_ht_r2ad_sig, shr_ht_r2ad_sig$log2FoldChange<0)
ELGvPCL_sig_up <- subset(ELGvPCL_sig, ELGvPCL_sig$log2FoldChange<0)
ELGvPCL_sig_down <- subset(ELGvPCL_sig, ELGvPCL_sig$log2FoldChange>0)
EGLht_up <- shr_ht_r2ad_up[rownames(shr_ht_r2ad_up) %in%rownames(ELGvPCL_sig_up), ]$gene
nrow(shr_ht_r2ad_up)
length(EGLht_up)
EGLht_down <- shr_ht_r2ad_down[rownames(shr_ht_r2ad_down) %in%rownames(ELGvPCL_sig_down), ]$gene
nrow(shr_ht_r2ad_down)
length(EGLht_down)
## 647 out of 1863 genes upregulated genes in Het organoids are also upregulated in EGL compared to PCL/Bulk
## 639 out of 1958 genes downregulated genes in Het organoids are also downregulated in EGL compared to PCL/Bulk

ht_unique_res_up <- subset(ht_unique_res, ht_unique_res$log2FoldChange>0)
ht_unique_res_down <- subset(ht_unique_res, ht_unique_res$log2FoldChange<0)

EGLuphtup_un <- ht_unique_res_up[rownames(ht_unique_res_up) %in%rownames(ELGvPCL_sig_up), ]$gene
EGLdwnhtdwn_un <- ht_unique_res_down[rownames(ht_unique_res_down) %in%rownames(ELGvPCL_sig_down), ]$gene
View(data.frame(EGLuphtup_un))
nrow(ht_unique_res)
length(EGLuphtup_un)
length(EGLdwnhtdwn_un)

## 319 out of 1782 unique DGE genes in Het organoids are also upregulated in EGL compared to PCL/Bulk!!

pheatmap(al_mat_sub1[rownames(al_mat_sub1) %in% EGLuphtup_un , ], 
         scale="row", annotation_col = colannoALD, show_rownames = T, 
         show_colnames = F)

## Now lets find highset upregulation
ht_unique_res_up2 <- ht_unique_res_up[rownames(ht_unique_res_up) %in%rownames(ELGvPCL_sig_up), ]
ht_up  <- shr_ht_r2ad_up[rownames(shr_ht_r2ad_up) %in%rownames(ELGvPCL_sig_up), ]

View(data.frame(ht_unique_res_up2))
View(data.frame(ht_up))

pheatmap(al_mat_sub1[rownames(al_mat_sub1) %in% ht_up[ht_up$log2FoldChange>1, ]$gene, ], 
         scale="row", annotation_col = colannoALD, show_rownames = T, 
         show_colnames = F)

EGLht_up2  <- ELGvPCL_sig_up[rownames(ELGvPCL_sig_up) %in%rownames(shr_ht_r2ad_up), ]
View(data.frame(EGLht_up2))
length(EGLht_up2[EGLht_up2$log2FoldChange < c(-3), ]$gene)
pheatmap(al_mat_sub1[rownames(al_mat_sub1) %in% EGLht_up2[EGLht_up2$log2FoldChange < c(-3.5), ]$gene, ], 
         scale="row", annotation_col = colannoALD, show_rownames = T, cluster_rows = F,
         show_colnames = F)

View(data.frame(ht_unique_res_up2))
View(data.frame(ht_up))

pheatmap(al_mat_sub1[rownames(al_mat_sub1) %in% ht_up[ht_up$log2FoldChange>1, ]$gene, ], 
         scale="row", annotation_col = colannoALD, show_rownames = T, 
         show_colnames = F)


genes <- shr_ht_r2ad_up[order(shr_ht_r2ad_up$padj), ] %>% .[1:50, ] %>% .$gene

pheatmap(al_mat_sub1[rownames(al_mat_sub1) %in% genes, ], 
         scale="row", annotation_col = colannoALD, show_rownames = T, 
         show_colnames = F)

## Genes upregulated in ptch1+/- but not in EGL:
htup_un_only <- outersect1(rownames(ht_unique_res_up), rownames(ELGvPCL_sig_up))
htup_un_only_res <- subset(shr_ht_r2ad, rownames(shr_ht_r2ad) %in% htup_un_only)
htup_un_only2 <- outersect1(rownames(htup_un_only_res), rownames(ELGvPCL_sig_down))
htup_un_only2_res <- subset(shr_ht_r2ad, rownames(shr_ht_r2ad) %in% htup_un_only2)
pheatmap(al_mat_sub1[rownames(al_mat_sub1) %in% htup_un_only2_res[htup_un_only2_res$log2FoldChange>1, ]$gene, ], 
         scale="none", annotation_col = colannoALD, show_rownames = T, 
         show_colnames = F)
View(data.frame(htup_un_only2_res))


### Aldinger stuff
#Load dataset
aldsamp <- read_xlsx("/Volumes/Partitie4/RNAseq_analysis/aldmve/Aldsamps.xlsx")
file.exists(aldsamp$files)
aldsamptxi <- tximeta(aldsamp)

# Generate count matrix for gene counts
gse_ald <- summarizeToGene(aldsamptxi)
colData(gse_ald)
aldds <- DESeqDataSet(gse_ald, ~ sex + Age + Line + PCL + RL_SVZ + RL_VZ + EGL)
aldds <- DESeq(aldds,
               fitType = "parametric",
               minReplicatesForReplace = Inf)

#determine which fit is most approriate:
parald <- estimateDispersions(aldds, fitType = "parametric")
residualsparald <- abs(log(mcols(parald)$dispGeneEst) - log(mcols(parald)$dispFit))
locald <- estimateDispersions(aldds, fitType = "local")
residualslocald <- abs(log(mcols(locald)$dispGeneEst) - log(mcols(locald)$dispFit))
summary(residualsparald) # median is 1.221
summary(residualslocald) # median is 1.006
# Maintain parametric fit!
resultsNames(aldds)

aldds <- DESeq(aldds,
               fitType = "local",
               minReplicatesForReplace = Inf)

alvsd <- vst(aldds, blind=FALSE, fitType = "local")

###
ald_egl <- lfcShrink(aldds, coef="EGL_Not_vs_EGL")
ald_egl$gene <- mcols(aldds)$gene_name
ald_egl_sig <- subset(ald_egl, ald_egl$padj<0.05)

EGL_up50 <- ald_egl_sig[order(ald_egl_sig$log2FoldChange), ] %>% .[1:50, ]
pheatmap(ald_mat_sub1[rownames(ald_mat_sub1) %in% EGL_up50$gene, ], scale="row", annotation_col = colannoALD2)

EGL_up <- subset(ald_egl_sig, ald_egl_sig$log2FoldChange<0)

summary(rownames(EGL_up) %in% rownames(ht_unique_res_up))

resultsNames(aldds)
ald_rlvz <- lfcShrink(aldds, coef="RL_VZ_RL_vz_vs_Not")
ald_rlvz$gene <- mcols(aldds)$gene_name
ald_rlvz_sig <- subset(ald_rlvz, ald_rlvz$padj<0.05)

RLVZ_up50 <- ald_rlvz_sig[order(ald_rlvz_sig$log2FoldChange), ] %>% .[1:50, ]
pheatmap(ald_mat_sub1[rownames(ald_mat_sub1) %in% RLVZ_up50$gene, ], scale="row", annotation_col = colannoALD2)

RLVZ_up <- subset(ald_rlvz_sig, ald_rlvz_sig$log2FoldChange>0)

summary(rownames(RLVZ_up) %in% rownames(ht_unique_res_up))

ald_rlsvz <- lfcShrink(aldds, coef="RL_SVZ_RL_svz_vs_Not")
ald_rlsvz$gene <- mcols(aldds)$gene_name
ald_rlsvz_sig <- subset(ald_rlsvz, ald_rlsvz$padj<0.05)

RLSVZ_up50 <- ald_rlsvz_sig[order(ald_rlsvz_sig$log2FoldChange), ] %>% .[1:50, ]
pheatmap(ald_mat_sub1[rownames(ald_mat_sub1) %in% RLSVZ_up50$gene, ], scale="row", annotation_col = colannoALD2)

RLSVZ_up <- subset(ald_rlsvz_sig, ald_rlsvz_sig$log2FoldChange>0)

summary(rownames(RLSVZ_up) %in% rownames(ht_unique_res_up))

ald_pcl <- lfcShrink(aldds, coef="PCL_PCL_vs_Not")
ald_pcl$gene <- mcols(aldds)$gene_name
ald_pcl_sig <- subset(ald_pcl, ald_pcl$padj<0.05)

PCL_up50 <- ald_pcl_sig[order(ald_pcl_sig$log2FoldChange), ] %>% .[1:50, ]
pheatmap(ald_mat_sub1[rownames(ald_mat_sub1) %in% PCL_up50$gene, ], scale="row", annotation_col = colannoALD2)

PCL_up <- subset(ald_pcl_sig, ald_pcl_sig$log2FoldChange>0)

summary(rownames(PCL_up) %in% rownames(ht_unique_res_up))

ht_un_EGL <- subset(ht_unique_res_up,  rownames(ht_unique_res_up) %in% rownames(EGL_up))
ht_un_RLVZ <- subset(ht_unique_res_up, rownames(ht_unique_res_up) %in% rownames(RLVZ_up))
ht_un_RLSVZ <- subset(ht_unique_res_up,rownames(ht_unique_res_up) %in% rownames(RLSVZ_up))
ht_un_PCL <- subset(ht_unique_res_up,  rownames(ht_unique_res_up) %in% rownames(PCL_up))

genes <- c(ht_un_EGL$gene,
               ht_un_RLVZ$gene,
               ht_un_RLSVZ$gene,
               ht_un_PCL$gene)
genes2 <- unique(genes)

pheatmap(ald_mat_sub1[rownames(ald_mat_sub1) %in% genes2, ], cluster_rows=T, 
         show_rownames = F, scale="row", annotation_col = colannoALD2,
         show_colnames = F)

nrow(subset(shr_ht_r2ad_sig,  rownames(shr_ht_r2ad_sig) %in% rownames(EGL_up))) # 957/3821
nrow(subset(shr_ht_r2ad_sig,  rownames(shr_ht_r2ad_sig)%in% rownames(RLVZ_up))) # 945/3821
nrow(subset(shr_ht_r2ad_sig,  rownames(shr_ht_r2ad_sig) %in% rownames(RLSVZ_up))) # 803/3821
nrow(subset(shr_ht_r2ad_sig,  rownames(shr_ht_r2ad_sig)%in% rownames(PCL_up))) # 645/3821

notEGL <- c(ht_un_RLVZ$gene,
            ht_un_RLSVZ$gene,
            ht_un_PCL$gene)
EGLonly <- outersect1(ht_un_EGL$gene, notEGL)

notRLVZ <- c(ht_un_EGL$gene,
          ht_un_RLSVZ$gene,
          ht_un_PCL$gene)
RLVZonly <- outersect1(ht_un_RLVZ$gene,notRLVZ)

notRLSVZ <- c(ht_un_EGL$gene,
          ht_un_RLVZ$gene,
          ht_un_PCL$gene)
RLSVZonly <- outersect1(ht_un_RLSVZ$gene, notRLSVZ)

notPCL <- c(ht_un_EGL$gene,
         ht_un_RLVZ$gene,
         ht_un_RLSVZ$gene)
PCLonly <- outersect1(ht_un_PCL$gene, notPCL)

subgroupgenes <- c(EGLonly, RLVZonly, RLSVZonly, PCLonly)

pheatmap(ald_mat_sub1[subgroupgenes, ], cluster_rows=F, 
         show_rownames = F, scale="row", annotation_col = colannoALD2,
         show_colnames = F)

notEGL2 <- c(RLVZ_up$gene,
             RLSVZ_up$gene,
             PCL_up$gene)
EGLonly2 <- outersect1(EGL_up$gene, notEGL2)

notRLVZ2 <- c(EGL_up$gene,
              RLSVZ_up$gene,
              PCL_up$gene)
RLVZonly2 <- outersect1(RLVZ_up$gene,notRLVZ2)

notRLSVZ2 <- c(EGL_up$gene,
               RLVZ_up$gene,
               PCL_up$gene)
RLSVZonly2 <- outersect1(RLSVZ_up$gene, notRLSVZ2)

notPCL2 <- c(EGL_up$gene,
             RLVZ_up$gene,
             RLSVZ_up$gene)
PCLonly2 <- outersect1(PCL_up$gene, notPCL2)

region <- c(EGLonly2, RLVZonly2, RLSVZonly2, PCLonly2)

pheatmap(ald_mat_sub1[region, ], cluster_rows=F,
         show_rownames = F, scale="row", annotation_col = colannoALD2,
         show_colnames = F)

rowanno <- data.frame(c(rep("EGL", length(EGLonly2)), rep("RL_VZ", length(RLVZonly2)), rep("RL_SVZ", length(RLSVZonly2)), rep("PCL", length(PCLonly2))))
rownames(rowanno) <- region
colnames(rowanno) <- "Region"
region2 <- subset(region, region %in% rownames(mat_r2ad))
rowanno2 <- subset(rowanno, rownames(rowanno) %in% region2)
pheatmap(mat_r2ad[region2, ], scale="none", 
         annotation_row = rowanno2, annotation_col=ca, cluster_rows = F,
         show_rownames = F)

topVarGenes2 <- head(order(rowVars(mat_r2ad), decreasing = TRUE), 1000)
mat_r2ad_2 <- mat_r2ad[topVarGenes2 , ]
region3 <- subset(region, region %in% rownames(mat_r2ad_2))
rowanno3 <- subset(rowanno, rownames(rowanno) %in% region3)
pheatmap(mat_r2ad_2[region3, ], scale="none", 
         annotation_row = rowanno3, annotation_col=ca, cluster_rows = F,
         show_rownames = F)

mat_r2ad_3 <- mat_r2ad[shr_ht_r2ad_sig$gene , ]
region4 <- subset(region, region %in% rownames(mat_r2ad_3))
rowanno4 <- subset(rowanno, rownames(rowanno) %in% region4)
pheatmap(mat_r2ad_3[region4, ], scale="row", 
         annotation_row = rowanno4, annotation_col=ca, cluster_rows = T,
         show_rownames = F)


##### Other design####
aldds2 <- DESeqDataSet(gse_ald, ~ sex + Age + Line + region + region:Age)
aldds2 <- DESeq(aldds2,
                fitType = "parametric",
                minReplicatesForReplace = Inf)

#determine which fit is most approriate:
parald2 <- estimateDispersions(aldds2, fitType = "parametric")
residualsparald2 <- abs(log(mcols(parald2)$dispGeneEst) - log(mcols(parald2)$dispFit))
locald2 <- estimateDispersions(aldds2, fitType = "local")
residualslocald2 <- abs(log(mcols(locald2)$dispGeneEst) - log(mcols(locald2)$dispFit))
summary(residualsparald2) # median is 1.227
summary(residualslocald2) # median is 1.029

resultsNames(aldds2)
aldds2$region <- relevel(aldds2$region , ref="Bulk")
design(aldds2) <- ~ Age + sex + Line + region
aldds2 <- DESeq(aldds2,
                fitType = "local",
                minReplicatesForReplace = Inf)

EGL_res <- lfcShrink(aldds2, coef="region_EGL_vs_Bulk")
PCL_res <- lfcShrink(aldds2, coef="region_PCL_vs_Bulk")
RLVZ_res <- lfcShrink(aldds2, coef="region_RL_vz_vs_Bulk")
RLSVZ_res <- lfcShrink(aldds2, coef="region_RL_svz_vs_Bulk")

EGL_res$gene <- mcols(aldds2)$gene_name
PCL_res$gene <- mcols(aldds2)$gene_name
RLVZ_res$gene <- mcols(aldds2)$gene_name
RLSVZ_res$gene <- mcols(aldds2)$gene_name

## VST expression matrix
alvsd2 <- vst(aldds2, blind=FALSE, fitType = "local")
ald_mat <- assay(alvsd)
rownames(ald_mat) <- mcols(alddsall)$gene_name
ald_mat <- data.frame(ald_mat)
colnames(ald_mat)
ald_mat_sub1 <- subset(ald_mat, rowVars(ald_mat)!=0)
colannoALD2 <- data.frame(colData(aldds)$region, colData(aldds)$age_pcw)
rownames(colannoALD2) <- colnames(ald_mat_sub1)
colnames(colannoALD2) <- c("Region", "PCW")

EGL_sig <- subset(EGL_res, EGL_res$padj<0.05)
nrow(EGL_sig)
EGL_genes <- rownames(EGL_sig)
PCL_sig <- subset(PCL_res, PCL_res$padj<0.05)
nrow(PCL_sig)
PCL_genes <- rownames(PCL_sig)
RLVZ_sig <- subset(RLVZ_res, RLVZ_res$padj<0.05)
nrow(RLVZ_sig)
RLVZ_genes <- rownames(RLVZ_sig)
RLSVZ_sig <- subset(RLSVZ_res, RLSVZ_res$padj<0.05)
nrow(RLSVZ_sig)
RLSVZ_genes <- rownames(RLSVZ_sig)

EGL_unique <- outersect1(EGL_genes, c(PCL_genes, RLVZ_genes, RLSVZ_genes))
length(EGL_unique)
PCL_unique <- outersect1(PCL_genes, c(EGL_genes, RLVZ_genes, RLSVZ_genes))
length(PCL_unique)
RLVZ_unique <- outersect1(RLVZ_genes, c(PCL_genes, EGL_genes, RLSVZ_genes))
length(RLVZ_unique)
RLSVZ_unique <- outersect1(RLSVZ_genes, c(PCL_genes, RLVZ_genes, EGL_genes))
length(EGL_unique)

EGL_uniq_sig <- EGL_sig[EGL_unique, ]
EGL_up <- EGL_uniq_sig[EGL_uniq_sig$log2FoldChange>0, ]
EGL_dwn <- EGL_uniq_sig[EGL_uniq_sig$log2FoldChange<0, ]
PCL_uniq_sig <- PCL_sig[PCL_unique, ]
PCL_up <- PCL_uniq_sig[PCL_uniq_sig$log2FoldChange>0, ]
PCL_dwn <- PCL_uniq_sig[PCL_uniq_sig$log2FoldChange<0, ]
RLVZ_uniq_sig <- RLVZ_sig[RLVZ_unique, ]
RLVZ_up <- RLVZ_uniq_sig[RLVZ_uniq_sig$log2FoldChange>0, ]
RLVZ_dwn <- RLVZ_uniq_sig[RLVZ_uniq_sig$log2FoldChange<0, ]
RLSVZ_uniq_sig <- RLSVZ_sig[RLSVZ_unique, ]
RLSVZ_up <- RLSVZ_uniq_sig[RLSVZ_uniq_sig$log2FoldChange>0, ]
RLSVZ_dwn <- RLSVZ_uniq_sig[RLSVZ_uniq_sig$log2FoldChange<0, ]

EGL_ht_up <- subset(shr_ht_r2ad_up, rownames(shr_ht_r2ad_up) %in% rownames(EGL_up))
nrow(EGL_ht_up)
EGL_ht_dwn <- subset(shr_ht_r2ad_down, rownames(shr_ht_r2ad_down) %in% rownames(EGL_dwn))
nrow(EGL_ht_dwn)
PCL_ht_up <- subset(shr_ht_r2ad_up, rownames(shr_ht_r2ad_up) %in% rownames(PCL_up))
nrow(PCL_ht_up)
PCL_ht_dwn <- subset(shr_ht_r2ad_down, rownames(shr_ht_r2ad_down) %in% rownames(PCL_dwn))
nrow(PCL_ht_dwn)
RLVZ_ht_up <- subset(shr_ht_r2ad_up, rownames(shr_ht_r2ad_up) %in% rownames(RLVZ_up))
nrow(RLVZ_ht_up)
RLVZ_ht_dwn <- subset(shr_ht_r2ad_down, rownames(shr_ht_r2ad_down) %in% rownames(RLVZ_dwn))
nrow(RLVZ_ht_dwn)
RLSVZ_ht_up <- subset(shr_ht_r2ad_up, rownames(shr_ht_r2ad_up) %in% rownames(RLSVZ_up))
nrow(RLSVZ_ht_up)
RLSVZ_ht_dwn <- subset(shr_ht_r2ad_down, rownames(shr_ht_r2ad_down) %in% rownames(RLSVZ_dwn))
nrow(RLSVZ_ht_dwn)

pheatmap(mat_sig[EGL_ht_up$gene, ])


EGL_sig$change <- ifelse(EGL_sig$log2FoldChange>1, "Up", ifelse(EGL_sig$log2FoldChange<c(-2), "Down", "ns"))
EGL_ra <- EGL_sig[, 6:7]
rownames(EGL_ra) <- EGL_ra$gene
EGL_ra <- subset(EGL_ra, rownames(EGL_ra) %in% rownames(mat_r2ad))
mat_sig <- subset(mat_r2ad, rownames(mat_r2ad) %in%shr_ht_r2ad_sig$gene)
EGL_ra2 <- subset(EGL_ra, rownames(EGL_ra) %in% rownames(mat_sig))
pheatmap(mat_sig[EGL_ra2[EGL_ra2$change=="Up", ]$gene, ], scale="row", annotation_col=ca, show_rownames = T, main="EGL up")
pheatmap(mat_sig[EGL_ra2[EGL_ra2$change=="Down", ]$gene, ], scale="row", annotation_col=ca, show_rownames = F, main="EGL down")

PCL_sig$change <- ifelse(PCL_sig$log2FoldChange>1, "Up", ifelse(PCL_sig$log2FoldChange<c(-2), "Down", "ns"))
PCL_ra <- PCL_sig[, 6:7]
rownames(PCL_ra) <- PCL_ra$gene
PCL_ra <- subset(PCL_ra, rownames(PCL_ra) %in% rownames(mat_r2ad))
mat_sig <- subset(mat_r2ad, rownames(mat_r2ad) %in%shr_ht_r2ad_sig$gene)
PCL_ra2 <- subset(PCL_ra, rownames(PCL_ra) %in% rownames(mat_sig))
pheatmap(mat_sig[PCL_ra2[PCL_ra2$change=="Up", ]$gene, ], scale="row", annotation_col=ca, show_rownames = T, main="PCL up")
pheatmap(mat_sig[PCL_ra2[PCL_ra2$change=="Down", ]$gene, ], scale="row", annotation_col=ca, show_rownames = F, main="PCL down")

RLVZ_sig$change <- ifelse(RLVZ_sig$log2FoldChange>1, "Up", ifelse(RLVZ_sig$log2FoldChange<c(-2), "Down", "ns"))
RLVZ_ra <- RLVZ_sig[, 6:7]
rownames(RLVZ_ra) <- RLVZ_ra$gene
RLVZ_ra <- subset(RLVZ_ra, rownames(RLVZ_ra) %in% rownames(mat_r2ad))
mat_sig <- subset(mat_r2ad, rownames(mat_r2ad) %in%shr_ht_r2ad_sig$gene)
RLVZ_ra2 <- subset(RLVZ_ra, rownames(RLVZ_ra) %in% rownames(mat_sig))
pheatmap(mat_sig[RLVZ_ra2[RLVZ_ra2$change=="Up", ]$gene, ], scale="row", annotation_col=ca, show_rownames = T, main="RLVZ up")
pheatmap(mat_sig[RLVZ_ra2[RLVZ_ra2$change=="Down", ]$gene, ], scale="row", annotation_col=ca, show_rownames = F, main="RLVZ down")

RLSVZ_sig$change <- ifelse(RLSVZ_sig$log2FoldChange>1, "Up", ifelse(RLSVZ_sig$log2FoldChange<c(-2), "Down", "ns"))
RLSVZ_ra <- RLSVZ_sig[, 6:7]
rownames(RLSVZ_ra) <- RLSVZ_ra$gene
RLSVZ_ra <- subset(RLSVZ_ra, rownames(RLSVZ_ra) %in% rownames(mat_r2ad))
mat_sig <- subset(mat_r2ad, rownames(mat_r2ad) %in%shr_ht_r2ad_sig$gene)
RLSVZ_ra2 <- subset(RLSVZ_ra, rownames(RLSVZ_ra) %in% rownames(mat_sig))
pheatmap(mat_sig[RLSVZ_ra2[RLSVZ_ra2$change=="Up", ]$gene, ], scale="row", annotation_col=ca, show_rownames = T, main="RLSVZ up")
pheatmap(mat_sig[RLSVZ_ra2[RLSVZ_ra2$change=="Down", ]$gene, ], scale="row", annotation_col=ca, show_rownames = F, main="RLSVZ down")

RL_der <- intersect(RLSVZ_ra2[RLSVZ_ra2$change=="Up", ]$gene, RLVZ_ra2[RLVZ_ra2$change=="Up", ]$gene)
length(RL_der)
EGLvRL_der <- intersect(EGL_ra2[EGL_ra2$change=="Up", ]$gene, RL_der)
length(EGLvRL_der)
RL_der_dwn <- intersect(RLSVZ_ra2[RLSVZ_ra2$change=="Down", ]$gene, RLVZ_ra2[RLVZ_ra2$change=="Down", ]$gene)
length(RL_der_dwn)
EGLvRL_der_dwn <- intersect(EGL_ra2[EGL_ra2$change=="Down", ]$gene, RL_der)
length(EGLvRL_der_dwn)
pheatmap(mat_r2ad[EGLvRL_der, ], scale="row", annotation_col = ca, main="Genes up in EGL+RL")
pheatmap(mat_r2ad[RL_der_dwn, ], scale="row", annotation_col = ca, main="Genes down in EGL+RL")
pheatmap(mat_r2ad[PCL_ra2[PCL_ra2$change=="Up", ]$gene, ], scale="row")

shr_ht_r2ad_sig %>% .[.$gene %in% EGLvRL_der, ] %>% data.frame(.) %>% View(.)
shr_ht_r2ad_sig %>% .[.$gene %in% RL_der_dwn, ] %>% data.frame(.) %>% View(.)
shr_ht_r2ad_sig %>% .[.$gene %in% EGL_ra2[EGL_ra2$change=="Down", ]$gene, ] %>% data.frame(.) %>% View(.)
shr_ht_r2ad_sig %>% .[.$gene %in% PCL_ra2[PCL_ra2$change=="Up", ]$gene, ] %>% data.frame(.) %>% View(.)

deregulated_vEGL <- shr_ht_r2ad_sig %>% .[.$gene %in% EGL_ra2[EGL_ra2$change=="Down", ]$gene, ] %>% subset(. , .$log2FoldChange>0.585) %>% .$gene
length(deregulated_vEGL)


### Integrate organoid data with aldinger data using Combat-Seq
#set wd
setwd("/Volumes/Partitie4/RNAseq_analysis_2/")

#Load dataset
test <- read_xlsx("/Volumes/Partitie4/RNAseq_analysis_2/aldvsmve.xlsx")
#file.exists(test$files)
#View(test)
#Transcript quantification import and generation of RangedSummarizedExperiment
# Ensembl - Homo sapiens - release 97 is used as transcriptome
se <- tximeta(test)
# assayNames(se)
# rowRanges(se)
# seqinfo(se)
# edb <- retrieveDb(se)
# class(edb)

# Generate count matrix for gene counts
gse <- summarizeToGene(se)

#Generate dataset which excludes G8 and F3:
gse_avm <- gse
dds_avm <- DESeqDataSet(gse_avm, ~ Age + region + batch2 + Genotype)
dds_avmua <- dds_avm
dds_avmua <- DESeq(dds_avm,
                  fitType = "parametric",
                  minReplicatesForReplace = Inf)

library(sva)
countmatrix <- counts(dds_avm)
batch <- colData(dds_avm)$batch2
group <- colData(dds_avm)$phenotype

adjusted_counts <- ComBat_seq(countmatrix, batch=batch, group=group)
coldata <- colData(dds_avm)
dds_avma <- DESeqDataSetFromMatrix(adjusted_counts,
                                  colData=coldata,
                                  design= ~ Age + region + batch2 + Genotype)
dds_avma$phenotype <- relevel(dds_avma$Genotype, ref="WT")

# data exploration ----
dds_avma <- DESeq(dds_avma,
                 fitType = "parametric",
                 minReplicatesForReplace = Inf)

## VSTs
vsd_avm <- vst(dds_avmua, blind=FALSE, fitType = "parametric")
vsd_avma <- vst(dds_avma, blind=FALSE, fitType = "parametric")

# PCA:
mat_avm <- assay(vsd_avm)
mat_avma <- assay(vsd_avma)
pca_avm <- prcomp(t(mat_avm))
pca_avma <- prcomp(t(mat_avma))

dims_avm <- data.frame(PC= factor(paste0("PC",1:65), levels=c(paste0("PC",1:65)),ordered=T),
                   var_explained=(pca_avm$sdev)^2/sum((pca_avm$sdev)^2))
dims_avma <- data.frame(PC= factor(paste0("PC",1:65), levels=c(paste0("PC",1:65)),ordered=T),
                       var_explained=(pca_avma$sdev)^2/sum((pca_avma$sdev)^2))

head(dims_avm)
head(dims_avma)

pca_avm_df <- data.frame(pca_avm$x)
pca_avma_df <- data.frame(pca_avma$x)
rownames(pca_avm_df) <- colnames(assay(vsd_avm))
rownames(pca_avma_df) <- colnames(assay(vsd_avma))
pca_avm_df <- cbind(pca_avm_df, vsd_avm$Genotype, vsd_avm$Tissue, vsd_avm$region, vsd_avm$Age)
pca_avma_df <- cbind(pca_avma_df, vsd_avma$Genotype, vsd_avma$Tissue, vsd_avma$region, vsd_avma$Age)
colnames(pca_avm_df)[c(66, 67, 68, 69)] <- c("Genotype", "Tissue", "Region", "Age")
colnames(pca_avma_df)[c(66, 67, 68, 69)] <- c("Genotype", "Tissue", "Region", "Age")
pcaplot_avm<- ggplot(pca_avm_df, aes(x=PC1, y=PC4, colour=Age, label=rownames(pca_avm_df))) + 
  geom_point(size=3) +
  geom_text_repel(max.overlaps = Inf) +
  theme_classic()+
  ggtitle("PCA")
pcaplot_avm
pcaplot_avma<- ggplot(pca_avma_df, aes(x=PC1, y=PC2, colour=Region, label=rownames(pca_avma_df))) + 
  geom_point(size=3) +
  geom_text_repel(max.overlaps = Inf) +
  theme_classic()+
  ggtitle("PCA")
pcaplot_avm
pcaplot_avma


## Sepp GSVA
CB_GSVA <- test <- read_xlsx("/Volumes/Partitie4/RNAseq_analysis_2/GSVA_CB_dv2.xlsx")
hh <- as.list(CB_GSVA, na.rm=TRUE)
CB_GSVA2 <- test <- read_xlsx("/Volumes/Partitie4/RNAseq_analysis_2/GSVA_CB_dv3.xlsx")
hh2 <- as.list(CB_GSVA2, na.omit=TRUE)

Sepp_gsva <- gsva(assay(vsd_r2adb), hh2, method = "gsva", kcdf = "Gaussian",min.sz = 15, max.sz = 500, mx.diff = TRUE)
Sepp_gsva <- Sepp_gsva[c(20, 6, 14, 15, 16, 23, 24, 5,4,3, 22, 12, 11, 10, 9, 7, 1, 17, 18, 2, 21, 8, 19, 13), ]
pheatmap(Sepp_gsva, annotation_col = ca, cluster_rows = F)

sepp_fit <- lmFit(Sepp_gsva, mod) %>%
eBayes()
decideTests(sepp_fit, p.value=0.05) %>%
summary()
sepp_sig <- topTable(sepp_fit, coef=2, n=Inf)
sepp_sig

pvals <- data.frame(sepp_fit$F.p.value)
rownames(pvals) <- rownames(Sepp_gsva)
sepp_sig$Significant <- ifelse(sepp_sig$adj.P.Val<0.05, "Yes", "No")
ca3 <- data.frame(ca[, -2])
rownames(ca3) <-rownames(ca)
colnames(ca3) <- "Genotype"
pheatmap(Sepp_gsva, annotation_col = ca, cluster_rows = T, show_colnames = F, cellwidth=20, cellheight=20, treeheight_col=4, treeheight_row=4, border_color="white", legend_labels = c("", ""), )

GCP_genes <- data.frame(hh$GCP)
GCP_genes$hh.GCP <- ifelse(GCP_genes$hh.GC=="NA", NA, GCP_genes$hh.GCP)
GCP_genes <- na.omit(GCP_genes)
nrow(GCP_genes)
GC_diff2_genes <- data.frame(hh$GC_diff_2)
GC_diff2_genes$hh.GC_diff_2 <- ifelse(GC_diff2_genes$hh.GC_diff_2=="NA", NA, GC_diff2_genes$hh.GC_diff_2)
GC_diff2_genes <- na.omit(GC_diff2_genes)
nrow(GC_diff2_genes)

summary(GCP_genes$hh.GCP %in% rownames(shr_ht_r2ad_sig)) # 53 true, 72 false
summary(GC_diff2_genes$hh.GC_diff_2 %in% rownames(shr_ht_r2ad_sig)) #14 true, 17 FALSE

GCP_genes$Type <- rep("GCP", nrow(GCP_genes))
colnames(GCP_genes) <- c("Gene", "Type")
GC_diff2_genes$Type <- rep("GC_diff_2", nrow(GC_diff2_genes))
colnames(GC_diff2_genes) <- c("Gene", "Type")
rownames(GCP_genes) <- GCP_genes$hh.GCP
rownames(GC_diff2_genes) <- GC_diff2_genes$hh.GC_diff_2
GC_GCP <- rbind(GCP_genes, GC_diff2_genes)
GC_GCP <- subset(GC_GCP, GC_GCP$Gene!='ENSG00000007372'& GC_GCP$Gene!='ENSG00000164438')
rownames(GC_GCP)<- GC_GCP$Gene
GC_GCPa <- merge(data.frame(shr_ht_r2ad_sig), GC_GCP, by=0)

View(GC_GCP2)
GC_GCP1 <-GC_GCPa[order(GC_GCPa$log2FoldChange), ] 
GC_GCP2 <-GC_GCPa[order(GC_GCPa$log2FoldChange, decreasing = T), ] 
GC_GCP3 <- GC_GCP1[1:15, ]
GC_GCP4 <- GC_GCP2[1:15, ]
GCPs5 <- rbind(GC_GCP3, GC_GCP4)
rownames(GCPs5) <- GCPs5$Row.names

ra <- GCPs5[, c(7, 9)]
ra2 <- data.frame(ra[, -1])
rownames(ra2) <- ra$gene


pheatmap(mat_het[GCPs5$gene, ], scale="row",
         annotation_col = ca, cluster_rows = T, show_colnames = F, cellwidth = 10, cellheight = 10,
         treeheight_col=4, treeheight_row=4, border_color="white", legend_labels = c("", ""),
         annotation_row=ra2)

GCPs <- c("PAX6", "ZIC1", "ZIC2", "PCNA", "MKI67", "DLGAP5", "CCND1", "NEUROD1", "NEUROD2", "CNTN1", "CNTN2", "UNC5C","TBR1", "EOMES", "CALB2", "GRM1", "DCX")

GCPs <- pheatmap(mat_het[GCPs, ], scale="row",
         annotation_col = ca, cluster_rows = T, cluster_cols = F,  show_colnames = F, cellwidth = 10, cellheight = 10,
         treeheight_col=4, treeheight_row=4, border_color="white", legend_labels = c("", ""))

SHH <- c("GLI1", "GLI2", "GLI3", "PTCH1")

pheatmap(mat_het[SHH, ], scale="row",
         annotation_col = ca, cluster_rows = T, cluster_cols = F, show_colnames = F, cellwidth = 12, cellheight = 12,
         treeheight_col=4, treeheight_row=4, border_color="white", legend_labels = c("", ""))

## Genes in GCP list:
pheatmap(assay(vsd_r2adb)[rownames(assay(vsd_r2adb)) %in% hh2[["GCP"]], ], scale="row" )

View(data.frame(gcp_genes))

pheatmap(mat_r2adb[gcp_genes$gene, ], scale="row", annotation_col=ca)

View(data.frame(gc_genes))
pheatmap(mat_r2adb[gc_genes$gene, ], scale="row", annotation_col=ca)

gcp_genes <- subset(shr_ht_r2ad_sig, rownames(shr_ht_r2ad_sig) %in% hh2[["GCP"]]) %>% subset(. , .$log2FoldChange>0)
progenitor_genes <- subset(shr_ht_r2ad_sig, rownames(shr_ht_r2ad_sig) %in% hh2[["progenitor"]]) %>% subset(. , .$log2FoldChange>0)
VZ_neuroblast_1_genes <- subset(shr_ht_r2ad_sig, rownames(shr_ht_r2ad_sig) %in% hh2[["VZ_neuroblast_1"]]) %>% subset(. , .$log2FoldChange>0)
glioblast_genes <- subset(shr_ht_r2ad_sig, rownames(shr_ht_r2ad_sig) %in% hh2[["glioblast"]]) %>% subset(. , .$log2FoldChange>0)
NTZ_neuroblast_1_genes <- subset(shr_ht_r2ad_sig, rownames(shr_ht_r2ad_sig) %in% hh2[["NTZ_neuroblast_1"]]) %>% subset(. , .$log2FoldChange>0)
GC_diff_2_genes <- subset(shr_ht_r2ad_sig, rownames(shr_ht_r2ad_sig) %in% hh2[["GC_diff_2"]]) %>% subset(. , .$log2FoldChange<0)

GCP <- pheatmap(mat_r2adb[gcp_genes$gene, ], scale="row", annotation_col=ca, main= "Sepp GCP genes")
progenitor <- pheatmap(mat_r2adb[progenitor_genes$gene, ], scale="row", annotation_col=ca, main="Sepp progenitor genes")
VZ_neuroblast_1 <- pheatmap(mat_r2adb[VZ_neuroblast_1_genes$gene, ], scale="row", annotation_col=ca, main="Sepp VZ_neuroblast_1 genes")
glioblast <- pheatmap(mat_r2adb[glioblast_genes$gene, ], scale="row", annotation_col=ca, main="Sepp glioblast genes")
NTZ_neuroblast_1 <- pheatmap(mat_r2adb[NTZ_neuroblast_1_genes$gene, ], scale="row", annotation_col=ca, main="Sepp NTZ_neuroblast_1 genes")
GC_diff_2 <- pheatmap(mat_r2adb[GC_diff_2_genes$gene, ], scale="row", annotation_col=ca, main="Sepp GC_diff_2 genes")

grid.arrange(GCP[[4]], progenitor[[4]], VZ_neuroblast_1[[4]], NTZ_neuroblast_1[[4]], glioblast[[4]], GC_diff_2[[4]], ncol=2 )

glioblast

tumor_uniq <- read_xlsx("/Users/max/Documents/tumorgenes.xlsx")
tumor_uniq_b <- subset(tumor_uniq, tumor_uniq$Gene %in% rownames(mat_r2adb))
rownames(tumor_uniq_b) <- tumor_uniq_b$Gene
SHH <- subset(tumor_uniq_b, tumor_uniq_b$TumClass=="MB_SHH")

pheatmap(mat_r2adb[SHH$Gene, ], annotation_col = ca, scale="row", cluster_rows = T, show_rownames = T)

# FDA approved drug targets??
FDA <- read_xlsx("/Volumes/Partitie4/RNAseq_analysis/FDAgenes.xlsx")
View(FDA)

FDA_sig_ht <- subset(FDA, FDA$Ensembl%in%rownames(shr_ht_r2ad_sig))
FDA_sig_ht <- subset(FDA_sig_ht, FDA_sig_ht$Gene%in%rownames(mat_r2adb))
View(FDA_sig_ht)
FDA_drugs <- c("ERBB4",  "YAP1", "HDAC1", "HDAC9", "HMMR", "TOP2A",  "IDH1", "SMO", "SQLE", "HMGCR", "PCSK9", "FDFT1")
pheatmap(mat_r2adb[FDA_drugs, ], scale="row", annotation_col=ca3, cellwidth=12, cellheight=12, treeheight_col=4, treeheight_row=4, border_color="white", show_colnames = F, legend=F)

hetcts <- counts(dds_r2ad)
CBTNcts <- counts(CBTNdds)
TVG <- head(order(rowVars(assay(CBTNvsd)), decreasing = TRUE), 600)
CBTNcts2 <- CBTNcts[TVG, ]
hetcts2 <- hetcts[rownames(hetcts)%in%rownames(CBTNcts2), ]
CBTNcts3 <- CBTNcts2[rownames(CBTNcts2)%in%rownames(hetcts2), ]
CBTN_hets <- cbind(hetcts2, CBTNcts3)
nrow(CBTN_hets)
# View(data.frame(colData(dds_r2ad)))
# View(data.frame(colData(CBTNdds)))

batch <- c(colData(dds_r2ad)$batch, rep(3, ncol(CBTNcts2)))
group <- c(rep("WT", 6), rep( "Het", 5), colData(CBTNdds)$harmonized_diagnosis)
# library(sva)
adjusted_counts <- ComBat_seq(CBTN_hets,  batch=batch, group=NULL)
colnames(adjusted_counts) <- c(colnames(vsd_r2ad), colnames(CBTNvsd))
coldata2 <- data.frame(c(rep("Organoid", 11), rep("Tumour", ncol(CBTNdds))))
rownames(coldata2) <- colnames(adjusted_counts)
coldata2$genotype <- group
colnames(coldata2)[1] <- "Tissue"

CBTN_het_CBS <- DESeqDataSetFromMatrix(adjusted_counts,
                                  colData=coldata2,
                                  design= ~ genotype)
CBTN_het_CBS$genotype <- relevel(CBTN_het_CBS$genotype, ref="WT")

CBTN_het_CBS <- DESeq(CBTN_het_CBS,
                 fitType = "parametric",
                 minReplicatesForReplace = Inf)
CBTN_het_CBS_vsd <- vst(CBTN_het_CBS, blind=FALSE, fitType = "parametric")

#pheatmap(assay(CBTN_het_CBS_vsd), scale="row", show_rownames = F, show_colnames = F, annotation_col=coldata2)

# PCA:
pca_cbtn_ht <- prcomp(t(assay(CBTN_het_CBS_vsd)))

dims <- data.frame(PC= factor(paste0("PC",1:140), levels=c(paste0("PC",1:140)),ordered=T),
                   var_explained=(pca_cbtn_ht$sdev)^2/sum((pca_cbtn_ht$sdev)^2))
head(dims, 4)

pca_cbtn_ht_df <- data.frame(pca_cbtn_ht$x)
rownames(pca_cbtn_ht_df) <- colnames(assay(CBTN_het_CBS_vsd))
pca_cbtn_ht_df2 <- cbind(pca_cbtn_ht_df, CBTN_het_CBS_vsd$genotype, CBTN_het_CBS_vsd$Tissue)
colnames(pca_cbtn_ht_df2)[141:142] <- c("Genotype", "Tissue")
# Define the number of colors you want
nb.cols <- 24
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
ggplot(pca_cbtn_ht_df2, aes(x=PC1, y=PC2, colour=Genotype, shape=Tissue)) + 
  geom_point(size=3) +
  theme_classic()+
  facet_zoom(x = Genotype %in%c("WT", "Medulloblastoma, SHH-activated")) +
  scale_color_manual(values = mycolors) +
  ggtitle("PCA")


# Assessing sample similarity
sampleDists <- dist(t(assay(CBTN_het_CBS_vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(CBTN_het_CBS_vsd$genotype, CBTN_het_CBS_vsd$Tissue, sep="-")
colnames(sampleDistMatrix) <- NULL
sampleDistplot1 <- pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)




# Generate count matrix for gene counts
gse <- summarizeToGene(se)

#Generate dataset which including only WT and Het:
gse_r1 <- gse[, -c(4:10, 16:23)]
dds_r1 <- DESeqDataSet(gse_r1, ~ genotype)
dds_r1$genotype <- relevel(dds_r1$genotype, ref="WT")

# data exploration ----
dds_r1 <- DESeq(dds_r1,
                  fitType = "parametric",
                  minReplicatesForReplace = Inf)
#Varient stabalizing transformation
vsd_r1 <- vst(dds_r1, blind=FALSE, fitType = "parametric")


pca_ua_r1 <- prcomp(t(assay(vsd_r1)))

dims <- data.frame(PC= factor(paste0("PC",1:8), levels=c(paste0("PC",1:8)),ordered=T),
                   var_explained=(pca_ua_r1$sdev)^2/sum((pca_ua_r1$sdev)^2))
head(dims, 4)
elbow<- dims[1:8,]
screeplot <- ggplot(elbow, aes(x=PC,y=var_explained, group=1))+
  geom_point(size=4)+
  geom_line()+
  labs(title="Scree plot, wt and ht only, unadjusted") +
  ylab("Variance explained")+
  theme_classic()
screeplot


pca_ua_r1df <- data.frame(pca_ua_r1$x)
rownames(pca_ua_r1df) <- colnames(assay(vsd_r1))
pca_ua_r1df <- cbind(pca_ua_r1df, vsd_r1$genotype)
colnames(pca_ua_r1df)[9] <- c("genotype")

pca_ua_r1df_wt <- pca_ua_r1df[pca_ua_r1df$genotype=="WT", ]
pca_ua_r1df_het <- pca_ua_r1df[pca_ua_r1df$genotype=="Het", ]
pca_ua_r1df_hom <- pca_ua_r1df[pca_ua_r1df$genotype=="Hom", ]


ggplot(pca_ua_r1df, aes(x=PC1, y=PC2, colour=genotype)) + 
  geom_point(size=3) +
  theme_classic()+   # change axis limits
  geom_mark_ellipse(expand = 0.04 ,aes(fill=genotype), alpha=0.05,)+
  coord_cartesian(clip = 'off') +
  ggtitle("PCA")
pcaplotr1


WriteXLS::WriteXLS(ht_unique_up, ExcelFileName="/Users/max/HTunique.xls")
WriteXLS::WriteXLS(data.frame(ht_unique_up), ExcelFileName="/Users/max/HTunique.xls")
opposites <- rownames(subset(shr_ht_r2ad_up, rownames(shr_ht_r2ad_up)%in%rownames(shr_hm_r2_sig[shr_hm_r2_sig$log2FoldChange<0, ])))
length(opposites)
WriteXLS::WriteXLS(data.frame(opposites), ExcelFileName="/Users/max/Opposites.xls")
head(universe, 5)
opposites_bp <- clusterProfiler::enrichGO(gene= opposites,
                                          universe      = universe,
                                          OrgDb         = "org.Hs.eg.db",
                                          keyType = "ENSEMBL",
                                          ont           = "BP",
                                          pAdjustMethod = "BH",
                                          pvalueCutoff  = 0.05,
                                          qvalueCutoff  = 0.05,
                                          readable      = TRUE)
barplot(opposites_bp)
dotplot(opposites_bp)
library(enrichplot)
dotplot(opposites_bp)
edox <- setReadable(opposites_bp, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=as.list(shr_ht_r2ad[, 2]))
library(DOSE)
edox <- setReadable(opposites_bp, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
edox <- setReadable(opposites_bp, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox)

opposites2 <- rownames(subset(shr_ht_r2ad_down, rownames(shr_ht_r2ad_down)%in%rownames(shr_hm_r2_sig[shr_hm_r2_sig$log2FoldChange>0, ])))
length(opposites2)
WriteXLS::WriteXLS(data.frame(opposites2), ExcelFileName="/Users/max/opposites2.xls")

opposites2_bp <- clusterProfiler::enrichGO(gene= opposites2,
                                          universe      = universe,
                                          OrgDb         = "org.Hs.eg.db",
                                          keyType = "ENSEMBL",
                                          ont           = "BP",
                                          pAdjustMethod = "BH",
                                          pvalueCutoff  = 0.05,
                                          qvalueCutoff  = 0.05,
                                          readable      = TRUE)
barplot(opposites2_bp)
dotplot(opposites2_bp)
library(enrichplot)
dotplot(opposites2_bp)
edox2 <- setReadable(opposites2_bp, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox2, foldChange=as.list(shr_ht_r2ad[, 2]))
library(DOSE)
edox <- setReadable(opposites2_bp, 'org.Hs.eg.db', 'ENTREZID')
p2 <- cnetplot(edox, foldChange=geneList)
edox <- setReadable(opposites2_bp, 'org.Hs.eg.db', 'ENTREZID')


opposites_cc <- clusterProfiler::enrichGO(gene= opposites,
                                          universe      = universe,
                                          OrgDb         = "org.Hs.eg.db",
                                          keyType = "ENSEMBL",
                                          ont           = "CC",
                                          pAdjustMethod = "BH",
                                          pvalueCutoff  = 0.05,
                                          qvalueCutoff  = 0.05,
                                          readable      = TRUE)
barplot(opposites_cc)
dotplot(opposites_cc)
edox_cc <- setReadable(opposites_cc, 'org.Hs.eg.db', 'ENTREZID')
p3 <- cnetplot(edox_cc, foldChange=as.list(shr_ht_r2ad[, 2]))

opposites2_cc <- clusterProfiler::enrichGO(gene= opposites2,
                                           universe      = universe,
                                           OrgDb         = "org.Hs.eg.db",
                                           keyType = "ENSEMBL",
                                           ont           = "CC",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05,
                                           readable      = TRUE)
barplot(opposites2_cc)
dotplot(opposites2_cc)
edox_cc2 <- setReadable(opposites2_cc, 'org.Hs.eg.db', 'ENTREZID')
p4 <- cnetplot(edox_cc2, foldChange=as.list(shr_ht_r2ad[, 2]))

opposites_MF <- clusterProfiler::enrichGO(gene= opposites,
                                          universe      = universe,
                                          OrgDb         = "org.Hs.eg.db",
                                          keyType = "ENSEMBL",
                                          ont           = "MF",
                                          pAdjustMethod = "BH",
                                          pvalueCutoff  = 0.05,
                                          qvalueCutoff  = 0.05,
                                          readable      = TRUE)
barplot(opposites_MF)
dotplot(opposites_MF)
edox_MF <- setReadable(opposites_MF, 'org.Hs.eg.db', 'ENTREZID')
p5 <- cnetplot(edox_MF, foldChange=as.list(shr_ht_r2ad[, 2]))

opposites2_MF <- clusterProfiler::enrichGO(gene= opposites2,
                                           universe      = universe,
                                           OrgDb         = "org.Hs.eg.db",
                                           keyType = "ENSEMBL",
                                           ont           = "MF",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05,
                                           readable      = TRUE)
barplot(opposites2_MF)
dotplot(opposites2_MF)
edox_MF2 <- setReadable(opposites2_MF, 'org.Hs.eg.db', 'ENTREZID')
p6 <- cnetplot(edox_MF2, foldChange=as.list(shr_ht_r2ad[, 2]))

edox3 <- pairwise_termsim(edox)
treeplot(edox3)
GL2 <- shr_ht_r2ad
rownames(GL2) <- GL2$gene
GL2 <- data.frame(GL2[, 2])
heatplot(edox, foldChange=GL2)

hu_bp <- clusterProfiler::enrichGO(gene= ht_unique_up,
                                           universe      = universe,
                                           OrgDb         = "org.Hs.eg.db",
                                           keyType = "ENSEMBL",
                                           ont           = "BP",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05,
                                           readable      = TRUE)
dotplot(hu_bp)
setReadable(hu_bp, 'org.Hs.eg.db', 'ENTREZID') %>%
  cnetplot(., foldChange=as.list(shr_ht_r2ad[, 2]))
setReadable(hu_bp, 'org.Hs.eg.db', 'ENTREZID') %>%
  pairwise_termsim(.) %>%
  treeplot(.)
upsetplot(hu_bp)

hu_CC <- clusterProfiler::enrichGO(gene= ht_unique_up,
                                   universe      = universe,
                                   OrgDb         = "org.Hs.eg.db",
                                   keyType = "ENSEMBL",
                                   ont           = "CC",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
dotplot(hu_CC)
setReadable(hu_CC, 'org.Hs.eg.db', 'ENTREZID') %>%
  cnetplot(., foldChange=as.list(shr_ht_r2ad[, 2]))
setReadable(hu_CC, 'org.Hs.eg.db', 'ENTREZID') %>%
  pairwise_termsim(.) %>%
  treeplot(.)
upsetplot(hu_CC)

hu_MF <- clusterProfiler::enrichGO(gene= ht_unique_up,
                                   universe      = universe,
                                   OrgDb         = "org.Hs.eg.db",
                                   keyType = "ENSEMBL",
                                   ont           = "MF",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
dotplot(hu_MF)
setReadable(hu_MF, 'org.Hs.eg.db', 'ENTREZID') %>%
  cnetplot(., foldChange=as.list(shr_ht_r2ad[, 2]))
setReadable(hu_MF, 'org.Hs.eg.db', 'ENTREZID') %>%
  pairwise_termsim(.) %>%
  treeplot(.)
upsetplot(hu_MF, main="test")

hu_bp2 <- clusterProfiler::enrichGO(gene= ht_unique_down,
                                   universe      = universe,
                                   OrgDb         = "org.Hs.eg.db",
                                   keyType = "ENSEMBL",
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
dotplot(hu_bp2)
setReadable(hu_bp2, 'org.Hs.eg.db', 'ENTREZID') %>%
  cnetplot(., foldChange=as.list(shr_ht_r2ad[, 2]))
setReadable(hu_bp2, 'org.Hs.eg.db', 'ENTREZID') %>%
  pairwise_termsim(.) %>%
  treeplot(.)
upsetplot(hu_bp2)

hu_CC2 <- clusterProfiler::enrichGO(gene= ht_unique_down,
                                   universe      = universe,
                                   OrgDb         = "org.Hs.eg.db",
                                   keyType = "ENSEMBL",
                                   ont           = "CC",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
dotplot(hu_CC2)
setReadable(hu_CC2, 'org.Hs.eg.db', 'ENTREZID') %>%
  cnetplot(., foldChange=as.list(shr_ht_r2ad[, 2]))
setReadable(hu_CC2, 'org.Hs.eg.db', 'ENTREZID') %>%
  pairwise_termsim(.) %>%
  treeplot(.)
upsetplot(hu_CC2)

hu_MF2 <- clusterProfiler::enrichGO(gene= ht_unique_down,
                                   universe      = universe,
                                   OrgDb         = "org.Hs.eg.db",
                                   keyType = "ENSEMBL",
                                   ont           = "MF",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.05,
                                   readable      = TRUE)
dotplot(hu_MF2)
setReadable(hu_MF2, 'org.Hs.eg.db', 'ENTREZID') %>%
  cnetplot(., foldChange=as.list(shr_ht_r2ad[, 2]))
setReadable(hu_MF2, 'org.Hs.eg.db', 'ENTREZID') %>%
  pairwise_termsim(.) %>%
  treeplot(.)
upsetplot(hu_MF2, main="test")

GO_opo_up <- clusterProfiler::enrichGO(gene= oposites[oposites$log2FoldChange.x>0.5, ]$rn,
                                  universe      = universe,
                                  OrgDb         = "org.Hs.eg.db",
                                  keyType = "ENSEMBL",
                                  ont           = "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.01,
                                  qvalueCutoff  = 0.05,
                                  readable      = TRUE)

GO_opo_down <- clusterProfiler::enrichGO(gene= oposites[oposites$log2FoldChange.x<c(-0.5), ]$rn,
                                       universe      = universe,
                                       OrgDb         = "org.Hs.eg.db",
                                       keyType = "ENSEMBL",
                                       ont           = "BP",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.01,
                                       qvalueCutoff  = 0.05,
                                       readable      = TRUE)

dotplot(GO_opo_up)
dotplot(GO_opo_down)

GO_ht_uni <- clusterProfiler::enrichGO(gene= rownames(ht_uniqe),
                                         universe      = universe,
                                         OrgDb         = "org.Hs.eg.db",
                                         keyType = "ENSEMBL",
                                         ont           = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff  = 0.01,
                                         qvalueCutoff  = 0.05,
                                         readable      = TRUE)


