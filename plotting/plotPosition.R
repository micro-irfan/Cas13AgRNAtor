library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
library(outliers)

setwd("C:/Irfan/Personal/Masters/BS6204/Project")

file <- "combined_output.raw.training.csv"
data <- read.table(file, header = TRUE, row.names = 1, sep=',')
selected_columns <- c("Score1", "Score2", "Score3", "Pos")
df <- data[, selected_columns]

df$Pos <- df$Pos+2
df$mean<-rowMeans(df[1:3], na.rm = TRUE)
df$Seq_id <- rownames(df)
df <- separate(df, Seq_id, into = c("Gene", "id"), sep = "-")

gene_lengths <- c(
  "CLUC" = 1733,
  "GLUC" = 608,
  "KRAS" = 5847,
  "MALAT1" = 8883,
  "PPIB" = 1059
)

# Convert to a data frame
gene_lengths_df <- as.data.frame(gene_lengths)
gene_lengths_df$Gene <- rownames(gene_lengths_df)

# Rename the columns if needed
colnames(gene_lengths_df) <- c("Length", "Gene")

df <- df %>%
  left_join(gene_lengths_df, by = "Gene")

df$PosPercent <- df$Pos / df$Length 

df_long <- df %>%
  pivot_longer(
    cols = starts_with("Score"),
    names_to = "Score_name",
    values_to = "Score_value"
  )

p <- ggplot(df_long, aes(x = PosPercent , y = Score_value, color = interaction(Gene))) +
  labs(x = "Position on Transcript", y = "Normalized Expression") + 
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  geom_point(alpha = 0.3) +
  geom_smooth(span = 0.025 , method = "loess", formula = y~x)

p

ggsave("Position.plot.png", plot = p, width = 10, height = 4)
