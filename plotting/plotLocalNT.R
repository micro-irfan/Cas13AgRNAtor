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

print_heatplot <- function(file, NAME) {
  nuc <- read.table(file, header = TRUE, row.names = 1, sep=',')
  colnames(nuc) <- sub("X", "", colnames(nuc))
  new_colnames <- c(paste0("-", 25:1), seq(LENGTH,1,-1), paste0("+", 1:25))
  colnames(nuc) <- new_colnames
  nuc$Window <- sub("row_", "", rownames(nuc))
  nuc_melt<-melt(nuc, id.vars = "Window", variable.name = "Position")
  nuc_melt$Window <- as.numeric(nuc_melt$Window)
  nuc_melt <- nuc_melt %>% arrange(as.numeric(Window))
  
  p<-ggplot(nuc_melt, aes(Position, Window)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient2(low = "darkblue", high = "red", mid = "white", guide = "colourbar", aesthetics = "fill",limits = scale_limits) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),  panel.grid = element_blank()) +
    theme_minimal() +
    scale_x_discrete(breaks = custom_x_labels, labels = custom_x_labels) +
    ggtitle(NAME) +
    xlab("Position In Target Site [p]")
  
  return(p)
} 

p <- list()

file <-"./RNA_context2/pearsonr.Nucleotide_context_A_edited.csv" 
p[['A']] <- print_heatplot(file, 'A Context') #+ scale_fill_gradient(limits = scale_limits)

file <-"./RNA_context2/pearsonr.Nucleotide_context_C_edited.csv"
p[['C']] <- print_heatplot(file, 'C Context') # + scale_fill_gradient(limits = scale_limits)

file <-"./RNA_context2/pearsonr.Nucleotide_context_G_edited.csv"
p[['G']] <- print_heatplot(file, 'G Context') #+ scale_fill_gradient(limits = scale_limits)

file <-"./RNA_context2/pearsonr.Nucleotide_context_T_edited.csv"
p[['U']] <- print_heatplot(file, 'U Context') #+ scale_fill_gradient(limits = scale_limits)

file <-"./RNA_context2/pearsonr.Nucleotide_context_AT_edited.csv"
p[['A|U']] <- print_heatplot(file, 'A|U Context') #+ scale_fill_gradient(limits = scale_limits)

file <-"./RNA_context2/pearsonr.Nucleotide_context_CG_edited.csv"
p[['C|G']] <- print_heatplot(file, 'C|G Context') #+ scale_fill_gradient(limits = scale_limits)


pdf("./figures/LocalNTContext.pdf",  width=20, height = 12, useDingbats = F)
grid.arrange(arrangeGrob(grobs= list(p[['A']], p[['C']], p[['G']], p[['U']], p[['A|U']], p[['C|G']] ),ncols=1))
dev.off()
