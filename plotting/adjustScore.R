library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
library(outliers)

setwd("C:/Irfan/Personal/Masters/BS6204/Project")

file <- "adjusted_raw_data.training.csv"
df <- read.table(file, header = TRUE, row.names = 1, sep=',')
df$Seq_id <- rownames(df)
df <- separate(df, Seq_id, into = c("Gene", "id"), sep = "-")

grouped_data <- df %>% group_by(Gene)
shapiro_test_results <- grouped_data %>%
  summarize(shapiro_p_value = shapiro.test(mean)$p.value)
print(shapiro_test_results)

result <- grouped_data %>% summarize(mean_expression = mean(mean))

print(result)

df1 <- df %>%
  left_join(result, by = "Gene") %>% # Join the grouped data to the original data
  mutate(adjustMean = mean - mean_expression) # Center Score1

df1$adjustMean <- -df1$adjustMean 

df1$quartile <- cut(df1$adjustMean, breaks = quantile(df1$adjustMean), labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)

df1$Seq_Id <- paste(df1$Gene, df1$id, sep = "-")
rownames(df1) <- df1$Seq_Id
df1 <- subset(df1, select = -c(gene, Seq_Id, id))
write.csv(df1, file = "MeanAdjustedScores.csv", row.names = TRUE)
