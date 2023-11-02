library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(stringr)
library(Biostrings)
library(ggpubr)

setwd("C:/Irfan/Personal/Masters/BS6204/Project")

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

file <- "BS6204_input.10WINDOW.csv"
data <- read.table(file, header = TRUE, row.names = 1, sep=',')
data$crRNA <- unlist(lapply(substring(data$Input, 11, 38), reverse_complement))


GetLetterFreq = function(G){
  x <- G
  df_row_names <- rownames(G) 
  G <- DNAStringSet(G$crRNA)
  
  AU.freq <- letterFrequency( G , letters = c("AT") ,as.prob = F)
  GC.freq <- letterFrequency( G , letters = c("GC") ,as.prob = F)
  AU.prob <- letterFrequency( G , letters = c("AT") ,as.prob = T)
  GC.prob <- letterFrequency( G , letters = c("GC") ,as.prob = T)
  A.prob <- letterFrequency( G , letters = c("A") ,as.prob = T)
  C.prob <- letterFrequency( G , letters = c("C") ,as.prob = T)
  G.prob <- letterFrequency( G , letters = c("G") ,as.prob = T)
  U.prob <- letterFrequency( G , letters = c("T") ,as.prob = T)
  A.freq <- letterFrequency( G , letters = c("A") ,as.prob = F)
  C.freq <- letterFrequency( G , letters = c("C") ,as.prob = F)
  G.freq <- letterFrequency( G , letters = c("G") ,as.prob = F)
  U.freq <- letterFrequency( G , letters = c("T") ,as.prob = F)
  diNucleotide.prob <- dinucleotideFrequency(G, as.prob = T)
  diNucleotide.freq <- dinucleotideFrequency(G, as.prob = F)
  
  out <- cbind.data.frame(A.freq,C.freq,G.freq,U.freq,A.prob,C.prob,G.prob,U.prob,GC.freq,AU.freq,GC.prob,AU.prob,diNucleotide.freq,diNucleotide.prob)
  rownames(out) <- df_row_names
  colnames(out) <- c("A","C","G","T","pA","pC","pG","pT","G|C", "A|T", "pG|pC", "pA|pT", "AA","AC","AG","AT","CA","CC","CG","CT","GA","GC",
                     "GG","GT","TA","TC","TG","TT","pAA","pAC","pAG","pAT","pCA","pCC","pCG","pCT","pGA","pGC","pGG","pGT","pTA","pTC","pTG","pTT")
  
  idx = match(rownames(out) , rownames(x))
  
  if (all(rownames(out) ==  rownames(x))  ){
    # out$Quartile <- x$Quartile[idx]
    out$guideScore <- x$guideScore[idx]
  }
  
  return(out)
}

GetNTfreqType = function(x){
  if (nchar(x) == 1){
    return("nucleotide")
  }
  else if (nchar(x) == 2){
    return("di-nucleotide")
  }
  else if ( grepl("\\|",x) ){
    return("nucleotides")
  }
  else{
    stop("Exiting! unexpected input")
  }
}


LollipopPlot = function( FREQ, NAME = ""){
  
  c = cor(FREQ[,grep("p",colnames(FREQ), invert = T)],use = "complete.obs")
  print (c)
  df = melt(c["guideScore",grep("guideScore", colnames(c) ,invert = T)])
  rownames(df) <- gsub("\\.","|",rownames(df))
  colnames(df) <- "pearson_coefficient"
  df$Nucleotide <- rownames(df)    
  df$type <- sapply(rownames(df), FUN = GetNTfreqType)
  df$type <- factor( df$type, levels = c("nucleotide","nucleotides","di-nucleotide"))
  df$Nucleotide <- gsub("T","U",df$Nucleotide)
  
  g1 = ggdotchart(df, x = "Nucleotide", y = "pearson_coefficient",
                  color = "type",                               # Color by groups
                  palette = c( "#E7B800", "#FC4E07","#00AFBB"), # Custom color palette
                  sorting = "descending",                       # Sort value in descending order
                  add = "segments",                             # Add segments from y = 0 to dots
                  add.params = list(color = "lightgray", size = 2), # Change segment color and size
                  #group = "type",                               # Order by groups
                  dot.size = 7,                                 # Large dot size
                  label = df$Nucleotide,                        # Add mpg values as dot labels
                  # label = round(df$pearson_coefficient,1),
                  font.label = list(color = "white", size = 9, vjust = 0.5),        # Adjust label parameters
                  ggtheme = theme_pubr(),                        # ggplot2 theme
                  ylim=c(-0.25,0.25),
                  title = NAME)
  
  return(g1)
  
}


letterfreq.all <- GetLetterFreq(data)
g1 = LollipopPlot(FREQ = letterfreq.all , NAME = "Perfect Match")
g1



