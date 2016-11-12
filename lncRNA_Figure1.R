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
                  
#making chr plot
library(ggplot2)
keeps <- c("id","V1", "V2", "V3", "V4")
lncRNA_keeps <- lncRNA_all[keeps]
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chr31", "chrUn", "chrX")
chrs_N <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "Un", "X")    
chromSizes<-read.table("equCab2.chrom.sizes.txt", header=T, col.names=c("chr","size"),stringsAsFactors=FALSE) 
m <- ggplot(data=lncRNA_keeps, aes(V1)) + geom_bar(aes(V1,group=id,colour=id),stat = "count") + scale_x_discrete(limits = chrs, labels = chrs_N) + ylab("lncRNA count") + xlab("Chr") +
  geom_line(data=subset(chromSizes,chr %in% chrs), aes(x=chr, y= size / 80000, group=1), colour="blue")

m

##Figure 1C
#merge lncRNA_keeps with expression data
setwd("~/Dropbox/Horse_Transcriptome/downloads")
overallExpression <- read.table("dataSummary", header=T, stringsAsFactors=F)
lncRNA_exp <- merge(overallExpression, lncRNA_keeps, by.x="transcriptName",by.y="V4" )
lncRNA_exp <-lncRNA_exp[ ,c("V1","V2","V3","transcriptName","calcTPM","id")]
#get exon information from unfiltered bed
setwd("~/Dropbox/Horse_Transcriptome/downloads/allTissues_BED")
unfiltered_bed <- read.table("unfiltered_Alltissues_Assembly.bed", header=F, stringsAsFactors=F)
lncRNA_exp_exons <- merge(lncRNA_exp, unfiltered_bed, by.x="transcriptName",by.y="V4" )
lncRNA_exp_exons <-lncRNA_exp_exons[ ,c("transcriptName","calcTPM","V10","id")]
#Figure 1C:cumulative expression and exons with categories
ggplot(subset(lncRNA_exp_exons, V10 %in% c(1:10)), aes(id,group=V10,fill=V10)) + 
  geom_bar(aes(weight=calcTPM),position="stack") + xlab("Input") + 
  ylab("cumulative expression (TPM)") + scale_fill_gradientn(colours=rainbow(10),
                                                             breaks=c(1,2,3,4,seq(5,10,2),10)) +
  guides(fill=guide_legend(title="exon number"))
