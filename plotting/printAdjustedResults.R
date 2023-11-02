library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)

setwd("C:/Irfan/Personal/Masters/BS6204/Project")

file <- "adjusted_raw_data.training.csv"
data <- read.table(file, header = TRUE, row.names = 1, sep=',')

kras <- data[data$gene == 'KRAS', ]
kras$mean2 <- -log2(kras$mean)
kras$quartile <- cut(kras$mean2, breaks = quantile(kras$mean2), labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)

Ymax=1.5
Ymin=-.5
f=0.2
b=500
LEN = 5724+28
g=ggplot(kras, aes(x=pos, y=mean2)) +
  geom_point(shape=20,  aes(color = quartile), size=2)  +
  geom_ribbon(aes(ymin=mean2-std, ymax=mean2+std), fill="blue2", alpha=0.2) +
  geom_smooth(span = 0.025 , method = "loess", formula = y~x, colour = "#525252" ) +
  theme_classic() +
  scale_color_manual(values = rev(c('#ca0020','#f4a582','#92c5de','#0571b0'))) +
  ylim(c(Ymin,Ymax)) + ylab("-log2(Normalized Kras Expression)") + xlab("guide match position [nt]") +
  # coord_fixed(ratio=(LEN/(Ymax-Ymin))*f) +
  scale_x_continuous(breaks = seq(0,LEN,b)) +
  theme(text = element_text(size=12)) +
  geom_hline(yintercept=0,linetype=3) 
g
ggsave("KRAS.adjusted.png", plot = g, width = 10, height = 4)
