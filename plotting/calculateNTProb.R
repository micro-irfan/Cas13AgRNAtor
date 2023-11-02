library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)

setwd("C:/Irfan/Personal/Masters/BS6204/Project")

file <- "BS6204_input.10WINDOW.csv"
data <- read.table(file, header = TRUE, row.names = 1, sep=',')

## Create Nucleotide Occurrences For Top20% Sequences

sorted_df <- data[order(-data$guideScore), ]
top_20_percent <- 0.20 * nrow(data)
top_20_percent_df <- head(sorted_df, n = top_20_percent)

LENGTH <- 28
WINDOW <- 10

position_range <- WINDOW + LENGTH + WINDOW
nucleotides <- c("A", "C", "G", "T")
sequences_20 <- top_20_percent_df$Input
count_matrix_20 <- matrix(0, nrow = position_range, ncol = length(nucleotides), dimnames = list(NULL, nucleotides))

for (seq_index in 1:length(sequences_20)) {
  sequence <- sequences_20[[seq_index]]
  for (position in 1:position_range) {
    nucleotide_at_position <- substr(sequence, position, position)
    count_matrix_20[position, nucleotide_at_position] <- count_matrix_20[position, nucleotide_at_position] + 1
  }
}

row_sums <- rowSums(count_matrix_20)
any(row_sums != length(sequences_20))

## Create Nucleotide Occurrences For All Sequences

sequences <- data$Input
count_matrix <- matrix(0, nrow = position_range, ncol = length(nucleotides), dimnames = list(NULL, nucleotides))

for (seq_index in 1:length(sequences)) {
  sequence <- sequences[[seq_index]]
  for (position in 1:position_range) {
    nucleotide_at_position <- substr(sequence, position, position)
    if (nucleotide_at_position == 'N'){
      next
    }
    count_matrix[position, nucleotide_at_position] <- count_matrix[position, nucleotide_at_position] + 1
  }
}

row_sums <- rowSums(count_matrix)
any(row_sums != length(sequences))

row_sums <- rowSums(count_matrix)
df_probabilities <- count_matrix / row_sums

## Calculate Delta Prob 
# this calculates a delta probability change (top - background)
# build fold change (=effect size)

row_sums_20 <- rowSums(count_matrix_20)
df_probabilities_20 <- count_matrix_20 / row_sums_20
delta_prob <- df_probabilities_20 - df_probabilities
effect_size <- log2(df_probabilities_20/df_probabilities)

out.favoured <- matrix(1,ncol = position_range, nrow = 4)
rownames(out.favoured ) = c("A","C","G","U")

out.disfavoured <- matrix(1,ncol = position_range, nrow = 4)
rownames(out.disfavoured ) = c("A","C","G","U")

for (i in 1:position_range) {
  observed_values <- count_matrix_20[i, ]
  expected_values <- count_matrix[i, ]
  expected_df <- df_probabilities[i, ]
  p.A = pbinom(observed_values['A'], size=sum(observed_values), prob = expected_df['A'], lower.tail = FALSE, log.p = F)
  p.C = pbinom(observed_values['C'], size=sum(observed_values), prob = expected_df['C'], lower.tail = FALSE, log.p = F)
  p.G = pbinom(observed_values['G'], size=sum(observed_values), prob = expected_df['G'], lower.tail = FALSE, log.p = F)
  p.U = pbinom(observed_values['T'], size=sum(observed_values), prob = expected_df['T'], lower.tail = FALSE, log.p = F)
  
  out.favoured["A",i] <-  p.A
  out.favoured["C",i] <-  p.C
  out.favoured["G",i] <-  p.G
  out.favoured["U",i] <-  p.U
  
  p.A = pbinom(observed_values['A'], size=sum(observed_values), prob = expected_df['A'], lower.tail = TRUE, log.p = F)
  p.C = pbinom(observed_values['C'], size=sum(observed_values), prob = expected_df['C'], lower.tail = TRUE, log.p = F)
  p.G = pbinom(observed_values['G'], size=sum(observed_values), prob = expected_df['G'], lower.tail = TRUE, log.p = F)
  p.U = pbinom(observed_values['T'], size=sum(observed_values), prob = expected_df['T'], lower.tail = TRUE, log.p = F)
  
  out.disfavoured["A",i] <-  p.A
  out.disfavoured["C",i] <-  p.C
  out.disfavoured["G",i] <-  p.G
  out.disfavoured["U",i] <-  p.U
}

## Adjust For Bonferroni
p.adjust( 0.004792735, method = "bonferroni", n = (28*4))

# any(out['G',] < 0.05)

transposed_delta_prob <- t(delta_prob)
rownames(transposed_delta_prob)[4]<- "U"

transposed_delta_prob <- transposed_delta_prob[,8:41]

colnames(transposed_delta_prob) <- c("-3","-2","-1",seq(LENGTH,1,-1),"+1","+2","+3")

xlbs <- colnames(transposed_delta_prob)
colnames(transposed_delta_prob) <- seq(1,length(xlbs),1) # assign surrogate column names not to confuse R
df = melt(transposed_delta_prob)
colnames(df) <- c("NT", "Pos","DeltaProbability") # rename columns
df$Pos <- rep(xlbs,  each = 4 ) # assign the original Columnnames to get the relative positions
df$Pos <- factor(df$Pos, levels = xlbs) # factorize the position
df$NT <- factor(df$NT, levels = c("A","C","G","U")) # factorize the NT

deltaProb_plot = ggplot( data = df, aes(x=Pos, y=DeltaProbability , fill = NT)) + geom_bar(stat="identity", position=position_dodge()) + theme_classic() + scale_fill_manual(values=c("#41ab5d", "#4292c6", "#fe9929","#cb181d")) + ylim(-.15,0.15) + ylab(expression(Delta~" Probability")) + xlab("guide position") +
  geom_vline(xintercept=seq(2,position_range,1)-0.5, linetype="solid", color = "#bdbdbd", size=0.1) + coord_fixed(position_range/2.4) + ggtitle( expression(Delta~"Probability - nt prevalence in top20% vs all guide")) +
  annotate( geom = "text" , y = -0.145, x = 1 , label = "disfavored" , angle = 90, hjust = "inward", color = "darkgrey", size = 2 ) +
  annotate( geom = "text" , y = 0.145, x = 1 , label = "favored" , angle = 90, hjust = "inward", color = "darkgrey", size = 2 )

p = out.disfavoured
transposed_effect_size <- t(effect_size)
rownames(transposed_effect_size)[4]<- "U"


transposed_delta_prob <- t(delta_prob)
p[transposed_delta_prob < 0] <- 1
p[transposed_delta_prob >= 0] <- 1
for ( i in 1:ncol(transposed_delta_prob)){
  for (k in 1:nrow(transposed_delta_prob)){
    if (transposed_delta_prob[k,i] > 0){
      p[k,i] <- -log10(out.disfavoured[k,i])
    }
    else{
      p[k,i] <- log10(out.disfavoured[k,i])
    }
  }
}

p <- p[,8:41]

colnames(p) <- c("-3","-2","-1",seq(LENGTH,1,-1),"+1","+2","+3")

xlbs <- colnames(p)
colnames(p) <- seq(1,length(xlbs),1)
df = melt(p)
colnames(df) <- c("NT", "Pos","pVal")
df$Pos <- rep(xlbs,  each = 4 )
df$Pos <- factor(df$Pos, levels = xlbs)
df$NT <- factor(df$NT, levels = c("A","C","G","U"))
pval_plot = ggplot( data = df, aes(x=Pos, y=pVal , fill = NT)) + geom_bar(stat="identity", position=position_dodge()) + theme_classic() + scale_fill_manual(values=c("#41ab5d", "#4292c6", "#fe9929","#cb181d")) + ylim(-4,4) + ylab("-log10(p)") + xlab("guide position") +
  geom_vline(xintercept=seq(2,position_range,1)-0.5, linetype="solid", color = "#bdbdbd", size=0.1) + coord_fixed(position_range/80) + ggtitle("Pval - nt prevalence in top20% v all guide") +
  annotate( geom = "text" , y = -4, x = 1 , label = "disfavored" , angle = 90, hjust = "inward", color = "darkgrey", size = 2 ) +
  annotate( geom = "text" , y = 4, x = 1 , label = "favored" , angle = 90, hjust = "inward", color = "darkgrey", size = 2 )

ggsave("pval_plot_nt.png", plot = pval_plot, width = 10, height = 4)


pdf("./figures/NtCorr.pdf", width = 7, height = 10, useDingbats = F)
grid.arrange(arrangeGrob(grobs = list( deltaProb_plot, pval_plot)   ,ncol=1))
dev.off()


#Reference
#p.U = pbinom(Top.Cnts["T"], size=sum(Top.Cnts), prob = NT.probs["T"], lower.tail = LowerTail, log.p = F)
