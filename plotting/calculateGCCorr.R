library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(stringr)

setwd("C:/Irfan/Personal/Masters/BS6204/Project")

file <- "BS6204_input.10WINDOW.csv"
data <- read.table(file, header = TRUE, row.names = 1, sep=',')

my_assert <- function(condition, error_message = "Assertion failed") {
  if (!condition) {
    stop(error_message)
  }
}

LENGTH = 28

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

# Define a function to calculate GC%
calculate_GC_percentage <- function(sequence) {
  sequence <- reverse_complement(substring(sequence, 11, 38))
  my_assert(nchar(sequence) == LENGTH, nchar(sequence))
  
  # Count the number of G and C bases
  count_G <- sum(str_count(sequence, "G"))
  count_C <- sum(str_count(sequence, "C"))
  
  # Calculate the GC%
  total_bases <- nchar(sequence)
  gc_percentage <- ((count_G + count_C) / total_bases)
  
  return(gc_percentage)
}

data$GC <- unlist(lapply(data$Input, calculate_GC_percentage))

fit <- lm(GC ~ guideScore, data = data)
r = summary(fit)$adj.r.squared
dat <- data[, c("guideScore", "GC")]
c <- cor(dat , use = "complete.obs" )[2]

summary_fit <- summary(fit)
p_value_guideScore <- summary_fit$coefficients["guideScore", "Pr(>|t|)"]

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


Y = sum(abs(range(data$guideScore , na.rm = T)))
X = max(data$GC , na.rm = T) - min(data$GC , na.rm = T) 
R = X/Y
My = max(data$guideScore , na.rm = T)
Mx = min(data$GC , na.rm = T)

g2 = ggplot(data, aes(x=GC, y=guideScore)) +
  geom_point(shape=20) +    # Use hollow circles
  geom_smooth(method=lm,   se=F, colour="#c82026")  +      
  geom_smooth(se = T, method = "loess", na.rm = T) +
  # geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = T) +
  theme_classic() +
  ylab("-log2(Normalized Expression)") + xlab("GC content") +
  annotate( geom = "text", label = paste0("rho = ",signif(c,2)), y = My-0.1, x = Mx, hjust = "inward") +
  annotate( geom = "text", label = paste0("R^2 = ",signif(r,2)), y = My-0.4, x = Mx, hjust = "inward") +
  annotate( geom = "text", label = paste0("p = ",signif(lmp(fit),2)), y = My-0.7, x = Mx, hjust = "inward") +  
  coord_fixed(ratio = R)
g2
