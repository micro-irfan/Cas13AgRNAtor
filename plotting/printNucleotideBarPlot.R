library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)

setwd("C:/Irfan/Personal/Masters/BS6204/Project")

### Plot Useless Plots From Casowary 
nt <- c("A", "T", "C", "G")

# Generate all possible combinations of length 2
combinations2 <- apply(expand.grid(nt, nt), 1, paste, collapse = "")
combinations3 <- apply(expand.grid(nt, nt, nt), 1, paste, collapse = "")

file <- "BS6204_input.10WINDOW.csv"
data <- read.table(file, header = TRUE, row.names = 1, sep=',')

LENGTH <- 28
WINDOW <- 10

reverse_complement <- function(dna_sequence) {
  # Define a mapping of nucleotides to their complements
  complement_mapping <- c(A = "T", T = "A", G = "C", C = "G")
  
  # Reverse the DNA sequence and replace nucleotides with their complements
  reversed_sequence <- rev(strsplit(dna_sequence, NULL)[[1]])
  complemented_sequence <- sapply(reversed_sequence, function(nucleotide) complement_mapping[nucleotide])
  
  # Combine the complemented nucleotides into a single string
  complemented_dna <- paste(complemented_sequence, collapse = "")
  
  return(complemented_dna)
}

data$crRNA <- unlist(lapply(substring(data$Input, 11, 38), reverse_complement))
sequences <- data$crRNA
count_matrix <- matrix(0, nrow = LENGTH, ncol = length(nt), dimnames = list(NULL, nt))

for (seq_index in 1:length(sequences)) {
  sequence <- sequences[[seq_index]]
  for (position in 1:LENGTH) {
    nucleotide_at_position <- substr(sequence, position, position)
    if (nucleotide_at_position == 'N'){
      next
    }
    count_matrix[position, nucleotide_at_position] <- count_matrix[position, nucleotide_at_position] + 1
  }
}

row_sums <- rowSums(count_matrix)
any(row_sums != length(sequences))

count_matrix_2 <- matrix(0, nrow = LENGTH-1, ncol = length(combinations2), dimnames = list(NULL, combinations2))

for (seq_index in 1:length(sequences)) {
  sequence <- sequences[[seq_index]]
  for (position in 1:(LENGTH-1)) {
    nucleotide_at_position <- substr(sequence, position, position+1)
    count_matrix_2[position, nucleotide_at_position] <- count_matrix_2[position, nucleotide_at_position] + 1
  }
}

row_sums_2 <- rowSums(count_matrix_2)
df_probabilities_2 <- count_matrix_2 / row_sums_2

data_to_plot <- df_probabilities_2[8,]  # Exclude the first column (Gene)
x_labels <- colnames(df_probabilities_2)  # Exclude the first column name (Gene)

plot_data_bad <- data.frame(
  Nucleotide = colnames(df_probabilities_2),  # Exclude the first column name (Gene)
  Percentage = df_probabilities_2[8,]  # Exclude the first column (Gene)
)


library(scales)
pbad <- ggplot(data = plot_data_bad, aes(x = reorder(Nucleotide, Percentage), y = Percentage, fill = Percentage)) +
  geom_bar(width=0.7, stat = "identity") +
  coord_flip() +
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grid lines
  scale_fill_gradient(low = "blue", high = "red") +
  xlab("Dinucleotide at crRNA Position 8") + 
  ylab("%")

ggsave("BarPlotPos8.png", plot = pbad, width = 6, height = 4)
pbad

result_dict <- list()
# Loop through sequences and guide scores
for (i in 1:length(sequences)) {
  dinucleotide <- substr(sequences[i], 8, 9)  # Extract the dinucleotide at position 8
  if (!is.null(result_dict[[dinucleotide]])) {
    result_dict[[dinucleotide]] <- c(result_dict[[dinucleotide]], data$guideScore[i])
  } else {
    result_dict[[dinucleotide]] <- list(data$guideScore[i])
  }
}

# Find the maximum number of guide scores
max_guide_scores <- max(sapply(result_dict, length))

# Pad data with NA for dinucleotides with fewer guide scores
for (key in names(result_dict)) {
  result_dict[[key]] <- c(result_dict[[key]], rep(NA, max_guide_scores - length(result_dict[[key]])))
}

# Convert the result_dict into a data frame
df <- data.frame(
  Dinucleotide = rep(names(result_dict), each = max_guide_scores),
  GuideScore = unlist(result_dict)
)

# Create a box plot
pbox <- ggplot(df, aes(x = Dinucleotide, y = GuideScore, fill = Dinucleotide)) +
  geom_boxplot(fill = "transparent") +
  theme_minimal() +
  ylab("Guide Score") +
  xlab('Dinucleotide at Position 8')
ggsave("BoxPlotPos8.png", plot = pbox, width = 6, height = 4)

pbox
