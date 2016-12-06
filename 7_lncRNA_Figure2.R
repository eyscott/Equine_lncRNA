##lncRNA figures
#have Table 1
#Figure 2
#C: chr plot
#combine all lncRNA db
setwd("~/Desktop/lncRNA")
lncRNA_all <- read.table("lncRNA_final_IDs", header=F, stringsAsFactors=F)
#making chr plot
library(ggplot2)
keeps <- c("V1", "V2", "V3", "V4","V5")
lncRNA_keeps <- lncRNA_all[keeps]
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chr31", "chrUn", "chrX")
chrs_N <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "Un", "X")    
chromSizes<-read.table("equCab2.chrom.sizes.txt", header=T, col.names=c("chr","size"),stringsAsFactors=FALSE) 
pdf("Fig1A.pdf")
ggplot(data=lncRNA_keeps, aes(V2)) + geom_bar(aes(V2,group=V1,fill=V1),stat = "count") + 
  scale_x_discrete(limits = chrs, labels = chrs_N) + ylab("lncRNA count") + xlab("Chr") +
  scale_fill_discrete(name  ="Input", 
                    guide = guide_legend(reverse = F),
                    labels=c("intergenic","known lncRNA","novel I","novel II","novel III")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  theme(axis.title = element_text(colour="black", size = 16)) +
  geom_line(data=subset(chromSizes,chr %in% chrs), aes(x=chr, y= size / 180000, group=1), colour="blue")

dev.off()

##Figure 2B
#merge lncRNA_keeps with expression data
setwd("~/Dropbox/Horse_Transcriptome/downloads")
overallExpression <- read.table("dataSummary", header=T, stringsAsFactors=F)
overallExpression$transcriptName=rownames(overallExpression)
lncRNA_exp <- merge(overallExpression, lncRNA_keeps, by.x="transcriptName",by.y="V5" )
lncRNA_exp <-lncRNA_exp[ ,c("V2","V3","V4","transcriptName","calcTPM","V1")]
#get exon information from unfiltered bed
unfiltered_bed <- read.table("allTissues_BED/unfiltered_Alltissues_Assembly.bed", header=F, stringsAsFactors=F)
lncRNA_exp_exons <- merge(lncRNA_exp, unfiltered_bed, by.x="transcriptName",by.y="V4" )
lncRNA_exp_exons <-lncRNA_exp_exons[ ,c("transcriptName","calcTPM","V10","V1.x")]
#Figure 2B:cumulative expression and exons with categories
library(RColorBrewer)
my.cols <- brewer.pal(6, "Set1")
pdf("Fig1C.pdf")
ggplot(subset(lncRNA_exp_exons, V10 %in% c(1:10)), aes(V1.x,group=V10,fill=V10)) + 
  geom_bar(aes(weight=calcTPM),position="stack") + xlab("Input") + 
  ylab("cumulative expression (TPM)") + scale_fill_gradientn(colours=my.cols,
                                                             breaks=c(1,2,3,4,5,6),
                                                          labels=c("1","2","3","4","5","6-10")) +
  scale_x_discrete(labels=c("intergenic","known lncRNA","novel I", "novel II", "novel III")) +
  guides(fill=guide_legend(title="exon number")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  theme(axis.title = element_text(colour="black", size = 16))
dev.off()

#Figure 2A: pie charts for overall RNAseq output for novel I,II,II, intergenic
#genes=those retained from hmmer and blastp
##Pie charts based on cumulative TPM
#get final lncRNA
setwd("~/Desktop/lncRNA")
all <- read.table("lncRNA_final_IDs", header=F, stringsAsFactors=F)
#attach expression values to these puppies
names(all)=c("id","V1","V2","V3")
#Attach coordinates to the expression values
overallExpression_coord <- merge(overallExpression,unfiltered_bed,by.x="transcriptName",by.y="V4")
overallExpression_coord <- overallExpression_coord[c("transcriptName","calcTPM","V1","V2","V3")]

all_exp <- merge(overallExpression_coord, all,by=c("V1","V2","V3"))  
keeps <- c("transcriptName","calcTPM","id")
all_exp <- all_exp[keeps]
novel_I_all_exp <- subset(all_exp, id %in% c("novel_I"))
novel_II_all_exp <- subset(all_exp, id %in% c("novel_II"))
novel_III_all_exp <- subset(all_exp, id %in% c("novel_III"))
intergenic_all_exp <- subset(all_exp, id %in% c("intergenic"))
known_all_exp <- subset(all_exp, id %in% c("known"))

#getting expression values for transcripts removed by filters
#Filter 1
F1_I <- read.table("F1_novel_I", header=F, stringsAsFactors=F)
F1_II <- read.table("F1_novel_II", header=F, stringsAsFactors=F)
F1_III <- read.table("F1_novel_III", header=F, stringsAsFactors=F)
F1_intergenic <- read.table("F1_intergenic", header=F, stringsAsFactors=F)
F1_known <- read.table("F1_known_ncRNA", header=F, stringsAsFactors=F)

#Filter 2
#F2_I <- read.table("F2_novel_I", header=F, stringsAsFactors=F)       ## no lines available in input
#F2_II <- read.table("F2_novel_II", header=F, stringsAsFactors=F)     ## no lines available in input
#F2_III <- read.table("F2_novel_III", header=F, stringsAsFactors=F)   ## no lines available in input
F2_intergenic <- read.table("F2_novel_intergenic", header=F, stringsAsFactors=F)
F2_known <- read.table("F2_known_ncRNA", header=F, stringsAsFactors=F)

#Filter 3
Pcoding <- read.table("P.bed", header=F, stringsAsFactors=F)
Pcoding_exp <- merge(Pcoding,overallExpression,by.x="V4",by.y="transcriptName")
F3_I <- subset(Pcoding_exp, id %in% c("novel_I"))
F3_I <- F3_I[ ,c("V4","length","calcTPM")]
F3_II <- subset(Pcoding_exp, id %in% c("novel_II"))
F3_II <- F3_II[ ,c("V4","length","calcTPM")]
F3_III <- subset(Pcoding_exp, id %in% c("novel_III"))
F3_III <- F3_III[ ,c("V4","length","calcTPM")]
F3_intergenic <- subset(Pcoding_exp, id %in% c("intergenic"))
F3_intergenic <- F3_intergenic[ ,c("V4","length","calcTPM")]
F3_known <- subset(Pcoding_exp, id %in% c("known"))
F3_known <- F3_known[ ,c("V4","length","calcTPM")]

#Filter 4
F4_I <- read.table("F4_novel_I", header=F, stringsAsFactors=F)
F4_II <- read.table("F4_novel_II", header=F, stringsAsFactors=F)
F4_III <- read.table("F4_novel_III", header=F, stringsAsFactors=F)
F4_intergenic <- read.table("F4_intergenic", header=F, stringsAsFactors=F)
F4_known <- read.table("F4_known", header=F, stringsAsFactors=F)

F4_I_exp <- merge(F4_I,overallExpression,by.x="V4",by.y="transcriptName")
F4_II_exp <- merge(F4_II,overallExpression,by.x="V4",by.y="transcriptName")
F4_III_exp <- merge(F4_III,overallExpression,by.x="V4",by.y="transcriptName")
F4_intergenic_exp <- merge(F4_intergenic,overallExpression,by.x="V4",by.y="transcriptName")
F4_known_exp <- merge(F4_known,overallExpression,by.x="V4",by.y="transcriptName")

F4_I_exp <- F4_I_exp[ ,c("V4","length","calcTPM")]
F4_II_exp <- F4_II_exp[ ,c("V4","length","calcTPM")]
F4_III_exp <- F4_III_exp[ ,c("V4","length","calcTPM")]
F4_intergenic_exp <- F4_intergenic_exp[ ,c("V4","length","calcTPM")]
F4_known_exp <- F4_known_exp[ ,c("V4","length","calcTPM")]

names(F3_I_exp)<-names(F1_I)
names(F3_II_exp)<-names(F1_II)
names(F3_III_exp)<-names(F1_III)
names(F3_intergenic_exp)<-names(F1_intergenic)

names(F4_I_exp)<-names(F1_I)
names(F4_II_exp)<-names(F1_II)
names(F4_III_exp)<-names(F1_III)
names(F4_intergenic_exp)<-names(F1_intergenic)

novel_I_rejects <- rbind(data.frame(id="novel_I_F1",F1_I),
                         #data.frame(id="novel_I_F2",F2_I),
                         data.frame(id="novel_I_F3",F3_I_exp),
                         data.frame(id="novel_I_F4",F4_I_exp))
novel_I_rejects_sub <- novel_I_rejects[ ,c("transcriptName","calcTPM","id")]

novel_II_rejects <- rbind(data.frame(id="novel_II_F1",F1_II),
                         #data.frame(id="novel_II_F2",F2_II),
                         data.frame(id="novel_II_F3",F3_II_exp),
                         data.frame(id="novel_II_F4",F4_II_exp))
novel_II_rejects_sub <- novel_II_rejects[ ,c("transcriptName","calcTPM","id")]

novel_III_rejects <- rbind(data.frame(id="novel_III_F1",F1_III),
                          #data.frame(id="novel_III_F2",F2_III),
                          data.frame(id="novel_III_F3",F3_III_exp),
                          data.frame(id="novel_III_F4",F4_III_exp))
novel_III_rejects_sub <- novel_III_rejects[ ,c("transcriptName","calcTPM","id")]

intergenic_rejects <- rbind(data.frame(id="intergenic_F1",F1_intergenic),
                           data.frame(id="intergenic_F2",F2_intergenic),
                           data.frame(id="intergenic_F3",F3_intergenic_exp),
                           data.frame(id="intergenic_F4",F4_intergenic_exp))
intergenic_rejects_sub <- intergenic_rejects[ ,c("transcriptName","calcTPM","id")]

known_rejects <- rbind(data.frame(id="known_F1",F1_known),
                           data.frame(id="known_F2",F2_known),
                           data.frame(id="known_F3",F3_known_exp),
                           data.frame(id="known_F4",F4_known_exp))
known_rejects_sub <- known_rejects[ ,c("transcriptName","calcTPM","id")]

#combine final lncRNA expression with rejects
novel_I_total_exp <-rbind(novel_I_all_exp,novel_I_rejects_sub)
novel_II_total_exp <-rbind(novel_II_all_exp,novel_II_rejects_sub)
novel_III_total_exp <-rbind(novel_III_all_exp,novel_III_rejects_sub)
intergenic_total_exp <-rbind(intergenic_all_exp,intergenic_rejects_sub)
known_total_exp <-rbind(known_all_exp,known_rejects_sub)

#Pie Chart based on cumulative TPM
library(RColorBrewer)
#my.cols <- brewer.pal(5, "Set3")
#my.cols# "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3"
my.cols <- c("#FFFFB3", "#BEBADA", "#80B1D3","#FB8072")

pdf("Fig2A_I.pdf")
ggplot(novel_I_total_exp,aes(x=factor(1),weight=calcTPM,fill=id)) + 
  geom_bar(width=1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_manual(values = my.cols,
                    labels=c("F1","F3","F4","lncRNA")) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
dev.off()

pdf("Fig2A_II.pdf")
ggplot(novel_II_total_exp,aes(x=factor(1),weight=calcTPM,fill=id)) + 
  geom_bar(width=1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_manual(values = my.cols,
                    labels=c("F1","F3","F4","lncRNA")) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
dev.off()

pdf("Fig2A_III.pdf")
ggplot(novel_III_total_exp,aes(x=factor(1),weight=calcTPM,fill=id)) + 
  geom_bar(width=1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_manual(values = my.cols,
                    labels=c("F1","F3","F4","lncRNA")) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
dev.off()

my.cols <- c("#FFFFB3","#8DD3C7", "#BEBADA", "#80B1D3","#FB8072")
pdf("Fig2A_intergenic.pdf")
ggplot(intergenic_total_exp,aes(x=factor(1),weight=calcTPM,fill=id)) + 
  geom_bar(width=1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_manual(values = my.cols,
                    labels=c("F1","F2","F3","F4","lncRNA")) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
dev.off()

pdf("Fig2A_known.pdf")
ggplot(known_total_exp,aes(x=factor(1),weight=calcTPM,fill=id)) + 
  geom_bar(width=1) + xlab(" ") + 
  ylab(" ") + coord_polar("y") + 
  guides(fill=guide_legend(title="composition")) +
  theme(legend.title = element_text(colour="black", size=18, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 16)) +
  theme(axis.text = element_text(colour="black", size = 14)) +
  scale_fill_manual(values = my.cols,
                    labels=c("F1","F2","F3","F4","lncRNA")) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())
dev.off()

#calculate some quick stats
lncRNA_exp <- merge(lncRNA_keeps, overallExpression, by.x="V5",by.y="transcriptName" )
lncRNA_exp <-lncRNA_exp[ ,c("V5","calcTPM","length","V1")]

novel_I_exp_l <- subset(lncRNA_exp, V1 %in% "novel_I")
novel_II_exp_l <- subset(lncRNA_exp, V1 %in% "novel_II")
novel_III_exp_l <- subset(lncRNA_exp, V1 %in% "novel_III")
intergenic_exp_l <- subset(lncRNA_exp, V1 %in% "intergenic")
known_exp_l <- subset(lncRNA_exp, V1 %in% "known")


mean_TPM_I<-mean(novel_I_exp_l[["calcTPM"]])
mean_TPM_II<-mean(novel_II_exp_l[["calcTPM"]])
mean_TPM_III<-mean(novel_III_exp_l[["calcTPM"]])
mean_TPM_intergenic<-mean(intergenic_exp_l[["calcTPM"]])
mean_TPM_known<-mean(known_exp_l[["calcTPM"]])

mean_length_I<-mean(novel_I_exp_l[["length"]])
mean_length_II<-mean(novel_II_exp_l[["length"]])
mean_length_III<-mean(novel_III_exp_l[["length"]])
mean_length_intergenic<-mean(intergenic_exp_l[["length"]])
mean_length_known<-mean(known_exp_l[["length"]])

total_length_I <-sum(novel_I_exp_l[["length"]])
total_length_II <-sum(novel_II_exp_l[["length"]])
total_length_III <-sum(novel_III_exp_l[["length"]])
total_length_intergenic <-sum(intergenic_exp_l[["length"]])
total_length_known <-sum(known_exp_l[["length"]])

sum_total <- sum(overallExpression[["length"]])
