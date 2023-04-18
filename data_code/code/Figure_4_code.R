library(readxl)
library(tidyverse)
library(rstatix)
library(DescTools)
library(ggplot2)
library(multcomp)
library(gridExtra)

df2 <- readRDS(file="Cyc_expr.rds")
WT <- df2 %>%
  subset(Line=="AH017")%>%
  ggplot(., aes(x=treatment, y=log2FC, fill=treatment))+
  stat_summary(fun.data="mean_sdl", 
               fun.args = list(mult=1), 
               geom="errorbar", color="grey15", alpha=0.5, width=0.3, size=0.5)+
  geom_bar(stat = "summary", fun.y = "mean", aes(fill=treatment), color="black", alpha=0.5, position=position_dodge()) + 
  geom_point(color="grey10", alpha=0.5, position=position_jitter(width=0.2)) +
  theme_classic()+
  scale_y_continuous( expand = expansion(mult = .2)) +
  scale_fill_brewer(palette="Set2")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  facet_wrap(~Gene, scale="free", ncol=2) +
  labs(x="")

HOM <- df2 %>%
  subset(Line=="H3")%>%
  ggplot(., aes(x=treatment, y=log2FC, fill=treatment))+
  stat_summary(fun.data="mean_sdl", 
               fun.args = list(mult=1), 
               geom="errorbar", color="grey15", alpha=0.5, width=0.3, size=0.5)+
  geom_bar(stat = "summary", fun.y = "mean", aes(fill=treatment), color="black", alpha=0.5, position=position_dodge()) + 
  geom_point(color="grey10", alpha=0.5, position=position_jitter(width=0.2)) +
  theme_classic()+
  scale_y_continuous( expand = expansion(mult = .2)) +
  scale_fill_brewer(palette="Set2")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  facet_wrap(~Gene, scale="free", ncol=2) +
  labs(x="")

WT + HOM

### DIV35 + 15 cease of treatment

df3 <- readRDS(file="Cyc_expr2_suppfile.rds")
df4 <- readRDS(file="Cyc_expr2.rds")
df4 <- df4*-1
rownames(df4) <- paste0(df3$Line, df3$treatment, 1:21)
df5 <- data.frame(df3[ , c(2, 3)])
rownames(df5) <- paste0(df3$Line, df3$treatment, 1:21)

#library(pheatmap)
pheatmap(df4, cluster_cols= F , cluster_rows = F ,scale="column", annotation_row = df5, show_colnames = T, show_rownames = T,
         cellwidth = 20, cellheight = 10, border_color = "white", gaps_row = 21)

cycpca <- prcomp(df4)

##
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

# PCA
ggbiplot(cycpca)+ 
  geom_point(size=3) +
  theme_classic()
