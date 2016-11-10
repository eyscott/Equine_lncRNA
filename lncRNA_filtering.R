##lncRNA pipeline R##
setwd("~/Desktop/lncRNA")
setwd("~/Dropbox/Horse_Transcriptome/downloads")
overallExpression <- read.table("dataSummary", header=T, stringsAsFactors=F)
overallExpression$transcriptName=rownames(overallExpression)

#load the different transcript categories
novel_I <- read.table("novelAnn/sup/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.sup", header=T, stringsAsFactors=F)
novel_II <- read.table("novelAnn/unsup.cons/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.unsup.cons", header=T, stringsAsFactors=F)
novel_III <- read.table("novelAnn/unsup.uncons.ORF/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.unsup.uncons.ORF", header=T, stringsAsFactors=F)
intergenic <- read.table("unsup", header=F, stringsAsFactors=F)
  
#attach expression values and gene lengths
novel_I <- merge(overallExpression, novel_I, by.x="transcriptName",by.y="transcript.ID" )
novel_II <- merge(overallExpression, novel_II, by.x="transcriptName",by.y="transcript.ID" )
novel_III <- merge(overallExpression, novel_III, by.x="transcriptName",by.y="transcript.ID" )
intergenic <- merge(overallExpression, intergenic, by.x="transcriptName",by.y="V1" )

#slim down to columns you want
keeps <- c("transcriptName","length", "calcTPM")
novel_I <- unique(novel_I[keeps])
novel_II <- unique(novel_II[keeps])
novel_III <- unique(novel_III[keeps])
intergenic <- unique(intergenic[keeps])

#filter out transcripts <200 nt and with TPM < 0.1
novel_I_f1 <- subset(novel_I, c(calcTPM > 0.1 & length > 200))
novel_II_f1 <- subset(novel_II, c(calcTPM > 0.1 & length > 200))
novel_III_f1 <- subset(novel_III, c(calcTPM > 0.1 & length > 200))
intergenic_f1 <- subset(intergenic, c(calcTPM > 0.1 & length > 200))
#Make each into a bed file to use with bedtools
setwd("~/Dropbox/Horse_Transcriptome/downloads/allTissues_BED")
unfiltered_bed <- read.table("unfiltered_Alltissues_Assembly.bed", header=F, stringsAsFactors=F)
novel_I_bed <- merge(novel_I_f1, unfiltered_bed, by.x="transcriptName",by.y="V4" )
novel_II_bed <- merge(novel_II_f1, unfiltered_bed, by.x="transcriptName",by.y="V4" )
novel_III_bed <- merge(novel_III_f1, unfiltered_bed, by.x="transcriptName",by.y="V4" )
intergenic_bed <- merge(intergenic_f1, unfiltered_bed, by.x="transcriptName",by.y="V4" )

#applying more stringent criteria to single exon transcripts
##V10 in bed is exon number
single_novel_I <- subset(novel_I_bed, c(calcTPM < 3 & length < 1000 & V10 <2))
single_novel_II <- subset(novel_II_bed, c(calcTPM < 3 & length < 1000 & V10 <2))
single_novel_III <- subset(novel_III_bed, c(calcTPM < 3 & length < 1000 & V10 <2))
single_intergenic <- subset(intergenic_bed, c(calcTPM < 3 & length < 1000 & V10 <2))
#remove this subset from bed file
require(dplyr)
novel_I_bed <- anti_join(novel_I_bed,single_novel_I, by="transcriptName")
novel_II_bed <- anti_join(novel_II_bed,single_novel_II, by="transcriptName")
novel_III_bed <- anti_join(novel_III_bed,single_novel_III, by="transcriptName")
intergenic_bed <- anti_join(intergenic_bed,single_intergenic, by="transcriptName")
#############

#calculate some quick stats
mean_TPM_I<-mean(novel_I_bed[["calcTPM"]])
mean_TPM_II<-mean(novel_II_bed[["calcTPM"]])
mean_TPM_III<-mean(novel_III_bed[["calcTPM"]])
mean_TPM_intergenic<-mean(intergenic_bed[["calcTPM"]])
mean_length_I<-mean(novel_I_bed[["length"]])
mean_length_II<-mean(novel_II_bed[["length"]])
mean_length_III<-mean(novel_III_bed[["length"]])
mean_length_intergenic<-mean(intergenic_bed[["length"]])
total_length_I <-sum(novel_I_bed[["length"]])
total_length_II <-sum(novel_II_bed[["length"]])
total_length_III <-sum(novel_III_bed[["length"]])
total_length_intergenic <-sum(intergenic_bed[["length"]])
sum_total <- sum(overallExpression[["length"]])
#calculating % of overlap
cov_novel_I <- (total_length_I/sum_total)
cov_novel_II <- (total_length_II/sum_total)
cov_novel_III <- (total_length_III/sum_total)
cov_intergenic<- (total_length_intergenic/sum_total)

#reorder for proper bed format
novel_I_bed <- novel_I_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_I_bed)[4]<-paste("V4")
novel_II_bed <- novel_II_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_II_bed)[4]<-paste("V4")
novel_III_bed <- novel_III_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_III_bed)[4]<-paste("V4")
intergenic_bed <- intergenic_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(intergenic_bed)[4]<-paste("V4")
#order chr properly
novel_I_bed <- novel_I_bed[with(novel_I_bed, order(V1, V2)), ]
novel_II_bed <- novel_II_bed[with(novel_II_bed, order(V1, V2)), ]
novel_III_bed <- novel_III_bed[with(novel_III_bed, order(V1, V2)), ]
intergenic_bed <- intergenic_bed[with(intergenic_bed, order(V1, V2)), ]
unfiltered_bed <- unfiltered_bed[with(unfiltered_bed, order(V1, V2)), ]
#write the bedfiles
setwd("~/Desktop/lncRNA")
write.table(novel_I_bed, "novel_I.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_bed, "novel_II.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_bed, "novel_III.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_bed, "intergenic.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(unfiltered_bed, "unfiltered.bed", row.names=F, col.names=F, quote=F, sep = "\t")
