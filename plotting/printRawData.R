library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)

setwd("C:/Irfan/Personal/Masters/BS6204/Project")

file1 <- "kras.scores.tsv"

kras<-read.delim(file1)

kras$mean<-rowMeans(kras[2:4], na.rm = TRUE)
kras$SD<-apply(kras[2:4], 1, sd, na.rm = TRUE)
kras$SE <- apply(kras[2:4], 1, function(row) sd(row, na.rm = TRUE) / sqrt(sum(!is.na(row))))
kras <- subset(kras, select = -X)

#goodlineplot

min_position<-0
max_position<-5724+28
breaks <- seq(0, max_position, by = 500)

pd <- position_dodge(0.1)
krasplot_lineplot_paper <- ggplot(kras, aes(x = Position, y = mean)) +
  geom_line(color = "blue", size=0.8) +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), size = 0.8, width = .1, color = "black") +
  labs(title = 'KRAS Line Plot With SE Error Bars', x = "Guide Position", y = NULL) +
  scale_y_continuous(name = "Normalized\nKRAS Expression") +  geom_line(position=pd) +
  geom_point(position=pd, size=1,alpha=0.5) +
  geom_hline(yintercept=1,linetype=3) +
  scale_x_continuous(limits = c(min_position, max_position), breaks = breaks) +
  theme_minimal()

krasplot_lineplot_paper

kras2_melt<-melt(kras, id.vars = c("Position","mean","SD"), measure.vars = c("Replicate1", "Replicate2","Replicate3"), variable.name = "Replicates")
krasplot_lineplot<-ggplot(kras2_melt, aes(x=Position, y=mean)) +
  geom_line(color="black", size=1) +
  geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD), fill="blue2", alpha=0.5) +
  labs(title = 'KRAS Line Plot With Geom_Ribbon', x = "Guide Position", y = NULL) +
  scale_y_continuous(name = "Normalized\nKRAS Expression") +  geom_hline(yintercept=1,linetype=3) +
  scale_x_continuous(limits = c(min_position, max_position), breaks = breaks) +
  theme_minimal()

krasplot_lineplot


krasplot_replicate<-ggplot(kras2_melt, aes(x=Position, y=value, group=Replicates)) +
  geom_line(aes(color=Replicates), show.legend = FALSE, size=.8, alpha=0.8) +
  labs(title = 'KRAS Line Plot With Replicate', x = "Guide Position", y = NULL) +
  scale_y_continuous(name = "Normalized\nKRAS Expression") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept=1,linetype=3) + 
  scale_x_continuous(limits = c(min_position, max_position), breaks = breaks) +
  theme_minimal() 

krasplot_replicate

p <- ggarrange(krasplot_lineplot, krasplot_lineplot_paper, krasplot_replicate,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)

ggsave("KRAS Raw Data.png", plot = p, width = 8, height = 8)

file2<-"malat1.scores.tsv"
malat<-read.delim(file2)

malat$mean<-rowMeans(malat[2:4], na.rm = TRUE)
malat$SD<-apply(malat[2:4], 1, sd, na.rm = TRUE)
malat$SE <- apply(malat[2:4], 1, function(row) sd(row, na.rm = TRUE) / sqrt(sum(!is.na(row))))

#goodlineplot

min_position<-0
max_position<-8652+28
breaks <- seq(0, max_position, by = 500)

pd <- position_dodge(0.1)
malatplot_lineplot_paper <- ggplot(malat, aes(x = Position, y = mean)) +
  geom_line(color = "blue", size=0.8) +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), size = 0.8, width = .1, color = "black") +
  labs(title = 'Malat1 Line Plot With SE Error Bars', x = "Guide Position", y = NULL) +
  scale_y_continuous(name = "Normalized\nMalat1 Expression") +  geom_line(position=pd) +
  geom_point(position=pd, size=1,alpha=0.5) +
  geom_hline(yintercept=1,linetype=3) +
  scale_x_continuous(limits = c(min_position, max_position), breaks = breaks) +
  theme_minimal()

malatplot_lineplot_paper

malat2_melt<-melt(malat, id.vars = c("Position","mean","SD"), measure.vars = c("Replicate1", "Replicate2","Replicate3"), variable.name = "Replicates")
malatplot_lineplot<-ggplot(malat2_melt, aes(x=Position, y=mean)) +
  geom_line(color="black", size=1) +
  geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD), fill="blue2", alpha=0.5) +
  labs(title = 'Malat1 Line Plot With Geom_Ribbon', x = "Guide Position", y = NULL) +
  scale_y_continuous(name = "Normalized\nMalat1 Expression") +  geom_hline(yintercept=1,linetype=3) +
  scale_x_continuous(limits = c(min_position, max_position), breaks = breaks) +
  theme_minimal()

malatplot_lineplot


malatplot_replicate<-ggplot(malat2_melt, aes(x=Position, y=value, group=Replicates)) +
  geom_line(aes(color=Replicates), show.legend = FALSE, size=.8, alpha=0.8) +
  labs(title = 'Malat1 Line Plot With Replicate', x = "Guide Position", y = NULL) +
  scale_y_continuous(name = "Normalized\nMalat1 Expression") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept=1,linetype=3) + 
  scale_x_continuous(limits = c(min_position, max_position), breaks = breaks) +
  theme_minimal() 

malatplot_replicate

p <- ggarrange(malatplot_lineplot, malatplot_lineplot_paper, malatplot_replicate,
          labels = c("A", "B", "C"),
          ncol = 1, nrow = 3)

ggsave("malat Raw Data.png", plot = p, width = 8, height = 8)
