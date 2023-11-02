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

selected_columns <- c("Score1", "Score2", "Score3")
df <- data[, selected_columns]

df$mean<-rowMeans(df[1:3], na.rm = TRUE)
df$meanAdjusted<--log2(df$mean)
df$SD<-apply(df[1:3], 1, sd, na.rm = TRUE)

df$Grubbs_Test_Result <- NA  
df$Shapiro_Result <- NA  


for (i in 1:nrow(df)) {
  row_data <- unlist(df[i, 1:3])  
  if (anyNA(row_data)) {
    next  # Skip this iteration and continue with the next one
  }
  outlier_test_result <- grubbs.test(row_data)
  df$Grubbs_Test_Result[i] <- outlier_test_result$p.value
  shapiro_results <-shapiro.test(row_data)
  df$Shapiro_Result[i] <-  shapiro_results$p.value
}

df$outlier <- df$Grubbs_Test_Result < 0.05

df$Seq_id <- rownames(df)
df <- separate(df, Seq_id, into = c("Gene", "id"), sep = "-")

grouped_data <- df %>% group_by(Gene)
result <- grouped_data %>%
  count(outlier)

result

df$Seq_id <- rownames(df)
df <- separate(df, Seq_id, into = c("Gene", "id"), sep = "-")

p <- ggplot(df, aes(x = mean, y = -log(Grubbs_Test_Result), color = Gene)) +
  geom_hline(yintercept=-log(0.05),linetype=3) +
  labs(x = "Normalized Expression", y = "-log10(pval)") + 
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1.25)) +
  geom_point()

p
ggsave("GrubbsTest.transformed.png", plot = p, width = 8, height = 8)

p1 <- ggplot(df, aes(x = Grubbs_Test_Result, y = mean, color = Gene)) +
  geom_vline(xintercept=0.05,linetype=3) +
  labs(x = "pval", y = "Normalized Expression") + 
  theme_minimal() +
  geom_point()

p1
ggsave("GrubbsTest.raw.png", plot = p1, width = 8, height = 8)


shapiro_test_results <- grouped_data %>%
  summarize(shapiro_p_value = shapiro.test(mean)$p.value)
print(shapiro_test_results)

result <- grouped_data %>% summarize(mean_expression = mean(mean))


df1 <- df %>%
  left_join(result, by = "Gene") %>% # Join the grouped data to the original data
  mutate(Score1_centered = Score1 - mean_expression,  # Center Score1
         Score2_centered = Score2 - mean_expression,  # Center Score2
         Score3_centered = Score3 - mean_expression)  # Center Score3
