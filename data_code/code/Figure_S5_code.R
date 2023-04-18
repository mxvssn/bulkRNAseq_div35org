library(readxl)
library(tidyverse)
library(rstatix)
library(DescTools)
library(ggplot2)
library(multcomp)
library(gridExtra)

df2 <- readRDS(file="iPSC_SHH_sig.rds")
df2 %>%
  ggplot(., aes(x=Line, y=FC))+
  stat_summary(fun.data="mean_sdl", 
               fun.args = list(mult=1), 
               geom="errorbar", color="grey15", width=0.3, size=0.5)+
  geom_bar(stat = "summary", aes(fill=Genotype), color="black", alpha=0.5, position=position_dodge()) + 
  geom_point(color="grey10", alpha=0.5, position=position_jitter(width=0.25), show.legend=FALSE) +
  theme_classic()+
  scale_fill_manual(values=c("#AED6F1", "#FFC300", "#FF5733"))+
  scale_y_continuous(trans= "log2", expand = expansion(mult = .2)) +
  scale_color_brewer(palette="Set2")+
  facet_wrap(~Gene, nrow=1) +
  labs(title="Gene expression in iPSCs", x="CRISPR Clones", y="Fold change")

df2 %>%
  subset(Gene=="GLI1") %>%
  drop_na(ddCT) %>%
  aov(ddCT ~ Line, data=.) %>% glht(., linfct = mcp(Line= "Dunnett")) %>% summary(.)

df2 %>%
  subset(Gene=="GLI1") %>%
  drop_na(ddCT) %>%
  aov(ddCT ~ Genotype, data=.) %>% glht(., linfct = mcp(Genotype= "Dunnett")) %>% summary(.)

## iPSC proliferation
edu <- readRDS(file="EdU_experiment.rds")
colnames(edu) <- c("Name", "Line", "Genotype", "Split", "G0_1", "G2", "S")
edu$Genotype <- factor(edu$Genotype, levels=c("WT", "Het", "Hom"), ordered=T)
edu_g <- gather(edu, key="Phase", value="Percentage", -Line, -Genotype, -Split, -Name)
edu_g$Percentage <- as.numeric(edu_g$Percentage)
edu_g$Phase <- factor(edu_g$Phase, levels=c("G0_1", "S","G2"), ordered=T)
edu_g$Line <- factor(edu_g$Line, levels=c("AH017", "A3", "B2", "B4", "C6", "A2", "B3", "H3"), ordered=T)
edu_g$Genotype <- factor(edu_g$Genotype, levels=c("WT", "Het", "Hom"), ordered=T)
edu_g$Split <- factor(edu_g$Split, levels=c("a", "b", "c", "d"), ordered=T)

edu_g %>% 
  subset(. , .$Phase%in%c("G0_1", "G2", "S")) %>%
  ggplot(.,  aes(x=Genotype, y=Percentage, group=Phase, color=Line)) +
  stat_summary(fun.data = mean_se, na.rm =F, position = "dodge", 
               geom = "bar", width = 0.8, size=1, color=c("#E41A1C", "#E41A1C", "#E41A1C", "#377EB8", "#377EB8", "#377EB8", "#4DAF4A", "#4DAF4A", "#4DAF4A"),
               fill=c("#E41A1C", "#E41A1C", "#E41A1C", "#377EB8", "#377EB8", "#377EB8", "#4DAF4A", "#4DAF4A", "#4DAF4A"), alpha=0.3) +
  geom_point(position=position_jitterdodge(dodge.width = 0.8), alpha=0.5)+
  stat_summary(fun.data="mean_se", 
               fun.args = list(mult=1), position=position_dodge(width=0.8),
               geom="errorbar", color="black", alpha=0.8, width=0.4, size=0.5)+
  theme_classic()+
  scale_fill_brewer(palette="Set1", labels = c("G0/1", "S", "G2")) +
  scale_color_brewer(palette="Set1") +
  labs(x="Clone", y="Percentage")