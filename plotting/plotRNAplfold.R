library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)

setwd("C:/Irfan/Personal/Masters/BS6204/Project")

LENGTH = 28
custom_x_labels <- c('-20', '-10', '-5', '-1', '+1', '+5', '+10', '+20')

scale_limits <- range(c(-0.25, 0.25))


file <-"pearsonr.RNAplfold.csv" 
nuc <- read.table(file, header = TRUE, row.names = 1, sep=',')
colnames(nuc) <- sub("X", "", colnames(nuc))
  
transposed_nuc <- t(nuc)

new_colnames <- c(paste0("-", 25:1), seq(LENGTH,1,-1), paste0("+", 1:25))
colnames(transposed_nuc) <- new_colnames

window_values <- rownames(transposed_nuc)
window_data <- data.frame(Window = window_values)
result_data <- cbind(window_data, transposed_nuc)

nuc_melt<-melt(result_data, id.vars = "Window", variable.name = "Position")
nuc_melt$Window <- as.numeric(nuc_melt$Window)
nuc_melt <- nuc_melt %>% arrange(as.numeric(Window))

my_colors <- c("red", "orange", "yellow", "lightblue", "darkblue")
my_colors <- c("darkblue", "lightblue",   "yellow", "red" )

p<-ggplot(nuc_melt, aes(Position, Window)) +
  geom_tile(aes(fill = value)) +
  #scale_fill_gradient2(low = "darkblue", high = "red", mid = "white", guide = "colourbar", aesthetics = "fill",limits = scale_limits) +
  scale_fill_gradientn(colors = my_colors, guide = "colourbar", aesthetics = "fill", limits = scale_limits) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  panel.grid = element_blank()) +
  theme_minimal() +
  scale_x_discrete(breaks = custom_x_labels, labels = custom_x_labels) +
  xlab("Position In Target Site [p]")
  
p
ggsave("RNAplfold.heatmap.png", plot = p, width = 10, height = 4)

p1<-ggplot(nuc_melt, aes(Position, Window)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "darkblue", high = "red", mid = "white", guide = "colourbar", aesthetics = "fill",limits = scale_limits) +
  #scale_fill_gradientn(colors = my_colors, guide = "colourbar", aesthetics = "fill", limits = scale_limits) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),  panel.grid = element_blank()) +
  theme_minimal() +
  scale_x_discrete(breaks = custom_x_labels, labels = custom_x_labels) +
  xlab("Position In Target Site [p]")

p1
ggsave("RNAplfold.heatmap.poor.png", plot = p1, width = 10, height = 4)
