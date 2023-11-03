library(ggplot2)
library(ggpubr)
library(readr)
library(gghighlight)
library(reshape2)
theme_update(plot.title = element_text(hjust = 0.5))

## Author: Denise Ng

#rna fold
file1 <- "RNAfold.results.outliers.csv"
rna<-read.csv(file1, sep = '\t')

rna$gq<-as.factor(rna$gq)
rna$dr<-as.factor(rna$dr)

#overplotting
overplot <- ggplot(rna, aes(x=mfe, y=guideScore, shape=dr, colour=gq, group=1)) + 
  geom_point(size=2) + 
  geom_smooth(method="lm") + 
  theme_bw() + 
  labs(x = "Minimum Free Energy (kcal/mol)", y="-log2 (Normalized Expression)")

#graph1
jitter<-position_jitter(width = 0.1, height = 0.1)
graph1 <- ggplot(rna, aes(x=mfe, y=guideScore, shape=dr, colour=gq, group=1)) + 
  geom_point(size=2,position=jitter) + 
  scale_shape_manual(values = c(1, 15)) + 
  scale_color_manual(values = c("blue", "red")) +
  geom_smooth(method="lm", color="black") +
  theme_bw() + 
  labs(x = "Minimum Free Energy (kcal/mol)", y="-log2 (Normalized Expression)")

#graph2
jitter<-position_jitter(width = 0.15, height = 0.15)
graph2 <- ggplot(rna, aes(x=mfe, y=guideScore, shape=dr, colour=gq, group=1)) + 
  geom_point(size=2,position=jitter) +
  scale_shape_manual(values = c(1, 15)) +
  scale_color_manual(values = c("#00AFBB", "red")) +
  geom_smooth(method="lm") + 
  theme_bw() + 
  labs(x = "Minimum Free Energy (kcal/mol)", y="-log2 (Normalized Expression)")

ggsave("RNAfold.overplot.png", plot = overplot, width = 10, height = 8)
ggsave("RNAfold.graph1.png", plot = graph1, width = 10, height = 8)
ggsave("RNAfold.graph2.png", plot = graph2, width = 10, height = 8)

melted<-melt(rna, id.vars = c("mfe", "guideScore"), measure.vars = c("gq", "dr"), variable.name = "Variable")


#boxplot

dr.labs<- c("Direct repeat = No", "Direct repeat = Yes")
names(dr.labs)<-c("0", "1")
gq.labs<- c("G-Quadruplex = No", "G-Quadruplex = Yes")
names(gq.labs)<-c("0", "1")

violin_mfe <- ggplot(rna, aes(x=mfe, y=guideScore, group=gq)) +
  geom_violin(aes(fill=gq),draw_quantiles = c(0.25, 0.5, 0.75)) +
  facet_grid(dr~ gq, labeller = labeller(dr=dr.labs, gq=gq.labs)) + 
  labs(x = "Minimum Free Energy (kcal/mol)", y="-log2 (Normalized Expression)")

box_mfe <- ggplot(rna, aes(x=mfe, y=guideScore)) + 
  geom_boxplot(fill="grey") + 
  facet_wrap(dr~gq, labeller = labeller(dr=dr.labs, gq=gq.labs)) +
  geom_jitter(size=1, position=position_jitter(0.2)) +
  labs(x = "Minimum Free Energy (kcal/mol)", y="-log2 (Normalized Expression)")

ggsave("RNAfold.violin.png", plot = violin_mfe, width = 10, height = 8)
ggsave("RNAfold.box.png", plot = box_mfe, width = 10, height = 8)


