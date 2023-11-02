library(ggplot2)
library(ggpubr)
library(readr)
library(gghighlight)
library(reshape2)

## Author: Denise Ng

file1 <- "corr.RNAhybrid.2.txt"
hybrid<-read.delim(file1, sep=",")

hybrid$start <- hybrid$start+1
hybrid$stop <- hybrid$stop+1
hybrid$Spearman <- hybrid$ScorrQ

#bbad heatmap
bad <- ggplot(data=hybrid, aes(x=start, y=stop, fill=Spearman)) + 
  geom_tile()+ggtitle("Spearman Correlation - Bad") +
  scale_y_reverse(position = "right",n.breaks=28,expand = c(0,0)) +
  scale_x_continuous(n.breaks=26,expand = c(0,0)) +
  theme_minimal() +
  labs(x = "Start Position [nt]", y="Stop Position [nt]")

#good heatmap
good <- ggplot(data=hybrid, aes(x=start, y=stop, fill=Spearman)) + 
  geom_tile() + 
  scale_fill_gradient2(low="darkblue", high="red2", mid="white", guide="colourbar", aesthetics="fill") +
  ggtitle("Spearman Correlation - Good") +
  scale_y_reverse(position = "right",n.breaks=28,expand = c(0,0)) +
  scale_x_continuous(n.breaks=26,expand = c(0,0))+labs(x = "Nucleotide position", y=" ") + 
  theme_minimal() +
  labs(x = "Start Position [nt]", y="Stop Position [nt]")

ggsave("RNAhybrid.good.png", plot = good, width = 10, height = 8)
ggsave("RNAhybrid.bad.png", plot = bad, width = 10, height = 8)


file2 <- "raw.RNAhybrid.2.csv"
raw<-read.csv(file2)

value<-max(abs(hybrid$SCorr))
highest<- hybrid[abs(hybrid$SCorr) == value,]

t_raw<-t(raw)
t_raw<-as.data.frame(t_raw)
colnames(t_raw)<-t_raw[1,]
t_raw<-t_raw[-1,]

file3<-"combined_output.removeOutliers.training.csv"
rna<-read.csv(file3)

new_df<-data.frame(
  column1=t_raw$`0-5 Mfe`,
  column2=rna$adjustedlog2FC,
  column3=rna$quartile
)

colnames(new_df)<-c("Mfe", "GuideScore", "Quartile")
new_df$Mfe <- as.numeric(new_df$Mfe)

raw <- ggplot(new_df, aes(x=Mfe, y=GuideScore, color = Quartile)) + 
  geom_point() + 
  geom_smooth(formula=y~x, method="lm", color='black') + 
  ggtitle("Guide scores for window (4-13)") + 
  labs(x = "Minimum Free Energy (kcal/mol)", y="-log2 (Normalized Expression)")

ggsave("RNAhybrid.raw.4.13.png", plot = raw, width = 10, height = 8)
