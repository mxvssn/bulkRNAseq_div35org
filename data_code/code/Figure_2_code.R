library(readxl)
library(tidyverse)
library(rstatix)
library(DescTools)
library(ggplot2)
library(multcomp)
library(gridExtra)
library(readxl)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(rstatix)

## Organodi growth
org_growth <- readRDS(file="organoid_growth.rds")
statistic <- org_growth %>% aov(ratio ~ genotype, data = .) %>% tukey_hsd()
ggplot(org_growth, aes(x=genotype, y=ratio)) + 
  geom_bar(stat="summary", fun = "mean", fill="azure", color="grey5", position="dodge")  + 
  stat_summary(fun.data="mean_sdl", 
               fun.args = list(mult=1), 
               geom="errorbar", color="black", width=0.5, size=1)+
  geom_point(aes(color=genotype), alpha=0.5, position=position_jitterdodge(
    jitter.width = 0.5,
    jitter.height = 0,
    dodge.width = 0.2,
    seed = NA
  ))  +
  stat_pvalue_manual(
    statistic,
    y.position = c(5.5, 6, 6.5),
    label = "p.adj.signif")+
  theme_classic()+
  scale_color_brewer(type="seq", palette="Set1") +
  labs(title= "Organoid growth", y="Growth", x="", color="Genotype") + theme(legend.position = "None")


## Gene expression at D35
df2 <- readRDS(file="D35_gene_expr.rds")

df2 %>%
  ggplot(., aes(x=Genotype, y=log2FC))+
  stat_summary(fun.data="mean_sdl", 
               fun.args = list(mult=1), 
               geom="errorbar", color="grey15", alpha=0.5, width=0.3, size=0.5)+
  geom_bar(stat = "summary", fun.y = "mean", aes(fill=Genotype), color="black", alpha=0.5, position=position_dodge()) + 
  geom_point(color="grey10", alpha=0.5, position=position_jitter(width=0.2)) +
  theme_classic()+
  scale_fill_manual(values=c("#AED6F1", "#FFC300", "#FF5733"))+
  scale_y_continuous( expand = expansion(mult = .2)) +
  scale_color_brewer(palette="Set2")+
  facet_wrap(~Gene, nrow=1) +
  labs(title="Gene expression at day 35", x="", caption="Relative expression to wildtype, 
                                                         normalized using ACTB and GAPDH")
