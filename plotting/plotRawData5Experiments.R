library(ggplot2)
library(ggpubr)
library(readr)
library(gghighlight)
library(reshape2)
theme_update(plot.title = element_text(hjust = 0.5))

# Author: Denise

#basic plots
cluc<-read.delim('C:/Users/Yong Quan/OneDrive/Documents/Denise/NTU/6203 Story telling/project/cluc.scores.tsv')
kras<-read.delim('C:/Users/Yong Quan/OneDrive/Documents/Denise/NTU/6203 Story telling/project/kras.scores.tsv')
malat<-read.delim('C:/Users/Yong Quan/OneDrive/Documents/Denise/NTU/6203 Story telling/project/malat1.scores.tsv')
ppib<-read.delim('C:/Users/Yong Quan/OneDrive/Documents/Denise/NTU/6203 Story telling/project/ppib.scores.tsv')
gluc<-read.delim('C:/Users/Yong Quan/OneDrive/Documents/Denise/NTU/6203 Story telling/project/gluc.scores.tsv')


#goodlineplot
cluc$mean<-rowMeans(cluc[2:4])
cluc$SD<-apply(cluc[2:4], 1, sd)
cluc2_melt<-melt(cluc, id.vars = c("Position","mean","SD"), measure.vars = c("Replicate1", "Replicate2","Replicate3"), variable.name = "Replicates")
clucplot2<-ggplot(cluc2_melt, aes(x=Position, y=mean))+geom_line(color="black", size=1)+geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD), fill="blue2", alpha=0.5)+theme_bw()+ggtitle("Cypridina luciferase")+geom_hline(yintercept=1,linetype=3)


#kras
kras$mean<-rowMeans(kras[2:4])
kras$SD<-apply(kras[2:4], 1, sd)
kras_melt<-melt(kras, id.vars = c("Position","mean","SD"), measure.vars = c("Replicate1", "Replicate2","Replicate3"), variable.name = "Replicates")
krasplot<-ggplot(kras_melt, aes(x=Position, y=mean))+geom_line(color="black")+geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD), fill="blue2", alpha=0.5)+theme_bw()+ggtitle("Kras")+geom_hline(yintercept=1,linetype=3)

#malat
malat$mean<-rowMeans(malat[2:4])
malat$SD<-apply(malat[2:4], 1, sd)
malat_melt<-melt(malat, id.vars = c("Position","mean","SD"), measure.vars = c("Replicate1", "Replicate2","Replicate3"), variable.name = "Replicates")
malatplot<-ggplot(malat_melt, aes(x=Position, y=mean))+geom_line(color="black")+geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD), fill="blue2", alpha=0.5)+theme_bw()+ggtitle("Malat1")+geom_hline(yintercept=1,linetype=3)

#ppib
ppib$mean<-rowMeans(ppib[2:4])
ppib$SD<-apply(ppib[2:4], 1, sd)
ppib_melt<-melt(ppib, id.vars = c("Position","mean","SD"), measure.vars = c("Replicate1", "Replicate2","Replicate3"), variable.name = "Replicates")
ppibplot<-ggplot(ppib_melt, aes(x=Position, y=mean))+geom_line(color="black")+geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD), fill="blue2", alpha=0.5)+theme_bw()+ggtitle("PPIB")+geom_hline(yintercept=1,linetype=3)

#gluc
gluc$mean<-rowMeans(gluc[2:4])
gluc$SD<-apply(gluc[2:4], 1, sd)
gluc_melt<-melt(gluc, id.vars = c("Position","mean","SD"), measure.vars = c("Replicate1", "Replicate2","Replicate3"), variable.name = "Replicates")
glucplot<-ggplot(gluc_melt, aes(x=Position, y=mean))+geom_line(color="black")+geom_ribbon(aes(ymin=mean-SD, ymax=mean+SD), fill="blue2", alpha=0.5)+theme_bw()+ggtitle("Gaussia luciferase")+geom_hline(yintercept=1,linetype=3)

ggarrange(clucplot2, krasplot, malatplot, ppibplot,glucplot,
          labels = c("A", "B", "C", "D", "E"),
          ncol = 2, nrow = 3)
   

#alternate lineplots with individual repliates
clucplot<-ggplot(cluc_melt, aes(x=Position, y=value, group=Replicates))+geom_line(aes(color=Replicates))+ggtitle("Cypridina luciferase") +theme(plot.title = element_text(hjust = 0.5))+geom_hline(yintercept=1,linetype=3)

bkrasplot<-ggplot(kras_melt, aes(x=Position, y=value, group=Replicates))+geom_line(aes(color=Replicates))+ggtitle("Kras") +theme(plot.title = element_text(hjust = 0.5))+geom_hline(yintercept=1,linetype=3)

bmalatplot<-ggplot(malat_melt, aes(x=Position, y=value, group=Replicates))+geom_line(aes(color=Replicates))+ggtitle("Malat") +theme(plot.title = element_text(hjust = 0.5))+geom_hline(yintercept=1,linetype=3)

bppibplot<-ggplot(ppib_melt, aes(x=Position, y=value, group=Replicates))+geom_line(aes(color=Replicates))+ggtitle("PPIB") +theme(plot.title = element_text(hjust = 0.5))+geom_hline(yintercept=1,linetype=3)

bglucbplot<-ggplot(gluc_melt, aes(x=Position, y=value, group=Replicates))+geom_line(aes(color=Replicates))+ggtitle("Gaussia luciferase") +theme(plot.title = element_text(hjust = 0.5))+geom_hline(yintercept=1,linetype=3)

ggarrange(clucplot, bkrasplot, bmalatplot, bppibplot,bglucbplot,
          labels = c("A", "B", "C", "D", "E"),
          ncol = 2, nrow = 3)

#nucleotide

nuc<-read.csv('C:/Users/Yong Quan/OneDrive/Documents/Denise/NTU/6203 Story telling/project/Nucleotide_context_A_edited.csv')
nuc_melt<-melt(nuc, id.vars = c("header"), variable.name = "Row")


#heatmap
ggplot(nuc_melt, aes(Row,header))+geom_tile(aes(fill=value))+scale_fill_gradient2(low="darkblue", high="red", mid="white",guide="colourbar", aesthetics="fill")+theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(nuc_melt, aes(Row,header))+geom_tile(aes(fill=value))+theme(axis.text.x = element_text(angle = 90, hjust = 1))


