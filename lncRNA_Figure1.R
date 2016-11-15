##lncRNA figures
#Will have Table 1
#Figure 1
#A) pipeline
#B: chr plot
#combine all lncRNA db
setwd("~/Desktop/lncRNA")
novel_I <- read.table("novel_I_f3.bed", header=F, stringsAsFactors=F)
novel_II <- read.table("novel_II_f3.bed", header=F, stringsAsFactors=F)
novel_III <- read.table("novel_III_f3.bed", header=F, stringsAsFactors=F)
intergenic <- read.table("intergenic_f3.bed", header=F, stringsAsFactors=F)
lncRNA_all <-rbind(data.frame(id="novel_I",novel_I),
                   data.frame(id="novel_II",novel_II),
                   data.frame(id="novel_III",novel_III),
                   data.frame(id="intergenic",intergenic))

lncRNA_all <- read.table("lncRNA_final8247_IDs", header=F, stringsAsFactors=F)
#making chr plot
library(ggplot2)
keeps <- c("V1", "V2", "V3", "V4","V5")
lncRNA_keeps <- lncRNA_all[keeps]
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chr31", "chrUn", "chrX")
chrs_N <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "Un", "X")    
chromSizes<-read.table("equCab2.chrom.sizes.txt", header=T, col.names=c("chr","size"),stringsAsFactors=FALSE) 
m <- ggplot(data=lncRNA_keeps, aes(V2)) + geom_bar(aes(V2,group=V1,fill=V1),stat = "count") + 
  scale_x_discrete(limits = chrs, labels = chrs_N) + ylab("lncRNA count") + xlab("Chr") +
  scale_fill_discrete(name  ="Input", 
                    guide = guide_legend(reverse = TRUE),
                    labels=c("intergenic","novel I","novel II","novel III")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  theme(axis.title = element_text(colour="black", size = 16)) +
  geom_line(data=subset(chromSizes,chr %in% chrs), aes(x=chr, y= size / 150000, group=1), colour="blue")

m

##Figure 1C
#merge lncRNA_keeps with expression data
setwd("~/Dropbox/Horse_Transcriptome/downloads")
overallExpression <- read.table("dataSummary", header=T, stringsAsFactors=F)
overallExpression$transcriptName=rownames(overallExpression)
lncRNA_exp <- merge(overallExpression, lncRNA_keeps, by.x="transcriptName",by.y="V5" )
lncRNA_exp <-lncRNA_exp[ ,c("V2","V3","V4","transcriptName","calcTPM","V1")]
#get exon information from unfiltered bed
setwd("~/Dropbox/Horse_Transcriptome/downloads/allTissues_BED")
unfiltered_bed <- read.table("unfiltered_Alltissues_Assembly.bed", header=F, stringsAsFactors=F)
lncRNA_exp_exons <- merge(lncRNA_exp, unfiltered_bed, by.x="transcriptName",by.y="V4" )
lncRNA_exp_exons <-lncRNA_exp_exons[ ,c("transcriptName","calcTPM","V10","V1.x")]
#Figure 1C:cumulative expression and exons with categories
ggplot(subset(lncRNA_exp_exons, V10 %in% c(1:10)), aes(V1.x,group=V10,fill=V10)) + 
  geom_bar(aes(weight=calcTPM),position="stack") + xlab("Input") + 
  ylab("cumulative expression (TPM)") + scale_fill_gradientn(colours=rainbow(10),
                                                             breaks=c(1,2,3,4,seq(5,10,2),10)) +
  scale_x_discrete(labels=c("intergenic","novel I", "novel II", "novel III")) +
  guides(fill=guide_legend(title="exon number")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  theme(axis.title = element_text(colour="black", size = 16))

#Figure 1D: pie charts for overall RNAseq output for novel I,II,II, intergenic
#genes=those retained from hmmer and blastp

setwd("~/Desktop/lncRNA")
all <- read.table("all_cats_PandnoP", header=F, stringsAsFactors=F)
novel_I_all <- subset(all, V1 %in% c("novel_I_lncRNA", "novel_I_genes"))
novel_II_all <- subset(all, V1 %in% c("novel_II_lncRNA", "novel_II_genes"))
novel_III_all <- subset(all, V1 %in% c("novel_III_lncRNA", "novel_III_genes"))
intergenic_all <- subset(all, V1 %in% c("intergenic_lncRNA", "intergenic_genes"))

##Pie charts based on counts
library(ggplot2)
ggplot(novel_I_all, aes(x=factor(1), fill=V1))+
  geom_bar(width = 1)+
  coord_polar("y") + ylab("") + xlab(" ") +
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_brewer(palette = "Set3",labels=c("genes","lncRNA"))
 
ggplot(novel_II_all, aes(x=factor(1), fill=V1))+
  geom_bar(width = 1)+
  coord_polar("y") + ylab("") + xlab(" ") +
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_brewer(palette = "Set3",labels=c("genes","lncRNA"))

ggplot(novel_III_all, aes(x=factor(1), fill=V1))+
  geom_bar(width = 1)+
  coord_polar("y") + ylab("") + xlab(" ") +
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_brewer(palette = "Set3",labels=c("genes","lncRNA"))

ggplot(intergenic_all, aes(x=factor(1), fill=V1))+
  geom_bar(width = 1)+
  coord_polar("y") + ylab("") + xlab(" ") +
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_brewer(palette = "Set3",labels=c("genes","lncRNA"))

##Pie charts based on cumulative TPM
setwd("~/Dropbox/Horse_Transcriptome/downloads")
overallExpression <- read.table("dataSummary", header=T, stringsAsFactors=F)
overallExpression$transcriptName=rownames(overallExpression)
all_exp <- merge(overallExpression, all, by.x="transcriptName",by.y="V5" )
keeps <- c("transcriptName","calcTPM","V1")
all_exp <- all_exp[keeps]
novel_I_all_exp <- subset(all_exp, V1 %in% c("novel_I_lncRNA", "novel_I_genes"))
novel_II_all_exp <- subset(all_exp, V1 %in% c("novel_II_lncRNA", "novel_II_genes"))
novel_III_all_exp <- subset(all_exp, V1 %in% c("novel_III_lncRNA", "novel_III_genes"))
intergenic_all_exp <- subset(all_exp, V1 %in% c("intergenic_lncRNA", "intergenic_genes"))

#Pie Chart based on cumulative TPM
ggplot(novel_I_all_exp,aes(x=factor(1),weight=calcTPM,fill=V1)) + 
  geom_bar(width=1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_brewer(palette = "Set3",labels=c("genes","lncRNA"))

ggplot(novel_II_all_exp,aes(x=factor(1),weight=calcTPM,fill=V1)) + 
  geom_bar(width=1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_brewer(palette = "Set3",labels=c("genes","lncRNA"))

ggplot(novel_III_all_exp,aes(x=factor(1),weight=calcTPM,fill=V1)) + 
  geom_bar(width=1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_brewer(palette = "Set3",labels=c("genes","lncRNA"))

ggplot(intergenic_all_exp,aes(x=factor(1),weight=calcTPM,fill=V1)) + 
  geom_bar(width=1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_brewer(palette = "Set3",labels=c("genes","lncRNA"))
