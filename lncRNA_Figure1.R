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

#binding all rejects from filters
F1_I <- read.table("F1_novel_I", header=F, stringsAsFactors=F)
F1_II <- read.table("F1_novel_II", header=F, stringsAsFactors=F)
F1_III <- read.table("F1_novel_III", header=F, stringsAsFactors=F)
F1_intergenic <- read.table("F1_novel_intergenic", header=F, stringsAsFactors=F)
F2_I <- read.table("F2_novel_I", header=F, stringsAsFactors=F)
F2_II <- read.table("F2_novel_II", header=F, stringsAsFactors=F)
F2_III <- read.table("F2_novel_III", header=F, stringsAsFactors=F)
F2_intergenic <- read.table("F2_novel_intergenic", header=F, stringsAsFactors=F)
F3_I <- read.table("F3_novel_I", header=F, stringsAsFactors=F)
F3_II <- read.table("F3_novel_II", header=F, stringsAsFactors=F)
F3_III <- read.table("F3_novel_III", header=F, stringsAsFactors=F)
F3_intergenic <- read.table("F3_intergenic", header=F, stringsAsFactors=F)
F3_I_ids <- data.frame(F3_I[ ,"V4"]) 
F3_II_ids <- data.frame(F3_II[ ,"V4"]) 
F3_III_ids <- data.frame(F3_III[ ,"V4"]) 
F3_intergenic_ids <- data.frame(F3_intergenic[ ,"V4"])
F3_I_exp <- merge(overallExpression, F3_I_ids, by.x="transcriptName",by.y="F3_I....V4.." )
F3_II_exp <- merge(overallExpression, F3_II_ids, by.x="transcriptName",by.y="F3_II....V4.." )
F3_III_exp <- merge(overallExpression, F3_III_ids, by.x="transcriptName",by.y="F3_III....V4.." )
F3_intergenic_exp <- merge(overallExpression, F3_intergenic_ids, by.x="transcriptName",by.y="F3_intergenic....V4.." )
F3_I_exp <- F3_I_exp[ ,c("transcriptName","length","calcTPM")]
F3_II_exp <- F3_II_exp[ ,c("transcriptName","length","calcTPM")]
F3_III_exp <- F3_III_exp[ ,c("transcriptName","length","calcTPM")]
F3_intergenic_exp <- F3_intergenic_exp[ ,c("transcriptName","length","calcTPM")]
names(F1_I)<-names(F3_I_exp)
names(F1_II)<-names(F3_II_exp)
names(F1_III)<-names(F3_III_exp)
names(F1_intergenic)<-names(F3_intergenic_exp)
names(F2_intergenic)<-names(F3_intergenic_exp)

novel_I_rejects <- rbind(data.frame(id="novel_I_F1",F1_I),
                         data.frame(id="novel_I_F3",F3_I_exp))
novel_I_rejects_sub <- novel_I_rejects[ ,c("transcriptName","calcTPM","id")]
names(novel_I_rejects_sub)[3]<-paste("V1")
novel_II_rejects <- rbind(data.frame(id="novel_II_F1",F1_II),
                         #data.frame(id="novel_II_F2",F2_II),
                         data.frame(id="novel_II_F3",F3_II_exp))
novel_II_rejects_sub <- novel_II_rejects[ ,c("transcriptName","calcTPM","id")]
names(novel_II_rejects_sub)[3]<-paste("V1")
novel_III_rejects <- rbind(data.frame(id="novel_III_F1",F1_III),
                          #data.frame(id="novel_III_F2",F2_III),
                          data.frame(id="novel_III_F3",F3_III_exp))
novel_III_rejects_sub <- novel_III_rejects[ ,c("transcriptName","calcTPM","id")]
names(novel_III_rejects_sub)[3]<-paste("V1")
intergenic_rejects <- rbind(data.frame(id="intergenic_F1",F1_intergenic),
                           data.frame(id="intergenic_F2",F2_intergenic),
                           data.frame(id="intergenic_F3",F3_intergenic_exp))
intergenic_rejects_sub <- intergenic_rejects[ ,c("transcriptName","calcTPM","id")]
names(intergenic_rejects_sub)[3]<-paste("V1")

novel_I_total_exp <-rbind(novel_I_all_exp,novel_I_rejects_sub)
novel_II_total_exp <-rbind(novel_II_all_exp,novel_II_rejects_sub)
novel_III_total_exp <-rbind(novel_III_all_exp,novel_III_rejects_sub)
intergenic_total_exp <-rbind(intergenic_all_exp,intergenic_rejects_sub)

#Pie Chart based on cumulative TPM
ggplot(novel_I_total_exp,aes(x=factor(1),weight=calcTPM,fill=V1)) + 
  geom_bar(width=1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_brewer(palette = "Set3",
                    labels=c("F1 rejects","F3 rejects","F4 rejects","lncRNA"))

ggplot(novel_II_total_exp,aes(x=factor(1),weight=calcTPM,fill=V1)) + 
  geom_bar(width=1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_brewer(palette = "Set3",
                    labels=c("F1 rejects","F3 rejects","F4 rejects","lncRNA"))

ggplot(novel_III_total_exp,aes(x=factor(1),weight=calcTPM,fill=V1)) + 
  geom_bar(width=1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_brewer(palette = "Set3",
                    labels=c("F1 rejects","F3 rejects","F4 rejects","lncRNA"))

ggplot(intergenic_total_exp,aes(x=factor(1),weight=calcTPM,fill=V1)) + 
  geom_bar(width=1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_brewer(palette = "Set3",
                    labels=c("F1 rejects","F2 reject","F3 rejects","F4 rejects","lncRNA"))

#calculate some quick stats
lncRNA_exp <- merge(lncRNA_keeps, overallExpression, by.x="V5",by.y="transcriptName" )
lncRNA_exp <-lncRNA_exp[ ,c("V5","calcTPM","length","V1")]

novel_I_exp_l <- subset(lncRNA_exp, V1 %in% "novel_I")
novel_II_exp_l <- subset(lncRNA_exp, V1 %in% "novel_II")
novel_III_exp_l <- subset(lncRNA_exp, V1 %in% "novel_III")
intergenic_exp_l <- subset(lncRNA_exp, V1 %in% "intergenic")


mean_TPM_I<-mean(novel_I_exp_l[["calcTPM"]])#6.08
mean_TPM_II<-mean(novel_II_exp_l[["calcTPM"]])#2.47
mean_TPM_III<-mean(novel_III_exp_l[["calcTPM"]])#1.35
mean_TPM_intergenic<-mean(intergenic_exp_l[["calcTPM"]])#1.08
mean_length_I<-mean(novel_I_exp_l[["length"]])#2192.6
mean_length_II<-mean(novel_II_exp_l[["length"]])#1956.1
mean_length_III<-mean(novel_III_exp_l[["length"]])#1933.9
mean_length_intergenic<-mean(intergenic_exp_l[["length"]])#1512.1
total_length_I <-sum(novel_I_exp_l[["length"]])#2014963
total_length_II <-sum(novel_II_exp_l[["length"]])#2153683
total_length_III <-sum(novel_III_exp_l[["length"]])#2417399
total_length_intergenic <-sum(intergenic_exp_l[["length"]])#7525619
sum_total <- sum(overallExpression[["length"]])#982439642
#calculating % of overlap
cov_novel_I <- (total_length_I/sum_total)#0.002
cov_novel_II <- (total_length_II/sum_total)#0.002
cov_novel_III <- (total_length_III/sum_total)#0.002
cov_intergenic<- (total_length_intergenic/sum_total)# 0.008
