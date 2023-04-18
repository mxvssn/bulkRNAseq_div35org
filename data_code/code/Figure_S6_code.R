library(readxl)
library(tidyverse)
library(rstatix)
library(DescTools)
library(ggplot2)
library(multcomp)
library(gridExtra)

setwd("/Users/max/Documents/PTCH1_manuscript/data_code/source")
df2 <- readRDS(file="D21_gene_expr.rds")
df2 %>%
  subset(Day==21)%>%
  ggplot(., aes(x=Genotype, y=FC))+
  stat_summary(fun.data="mean_sdl", 
               fun.args = list(mult=1), 
               position=position_dodge(width=0.85),
               geom="errorbar", color="grey15", alpha=0.5, width=0.3, size=0.5)+
  geom_bar(stat = "summary", fun.y = "mean", aes(fill=Genotype), color="black", alpha=0.5, position=position_dodge()) + 
  geom_point(aes(color=Genotype),alpha=0.5, position=position_jitterdodge(dodge.width=0.85, jitter.width = 0.5)) +
  theme_classic()+
  scale_fill_manual(values=c("#AED6F1", "#FFC300", "#FF5733"))+
  scale_y_continuous(trans= "log2", expand = expansion(mult = .2)) +
  scale_color_brewer(palette="Set2")+
  facet_wrap(~Gene, nrow=1) +
  labs(title="", x="")

#EN1
df2 %>%
  drop_na(avg_dCT) %>%
  .[.$Gene=="EN1", ]%>%
  aov(avg_dCT ~ Genotype, data=.) %>% glht(., linfct = mcp(Genotype= "Dunnett")) %>% summary(.)

#GBX2
df2 %>%
  drop_na(avg_dCT) %>%
  .[.$Gene=="GBX2", ]%>%
  aov(avg_dCT ~ Genotype, data=.) %>% glht(., linfct = mcp(Genotype= "Dunnett")) %>% summary(.)

#OTX2
df2 %>%
  drop_na(avg_dCT) %>%
  .[.$Gene=="OTX2", ]%>%
  aov(avg_dCT ~ Genotype, data=.) %>% glht(., linfct = mcp(Genotype= "Dunnett")) %>% summary(.)