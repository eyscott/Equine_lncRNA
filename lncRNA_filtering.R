##lncRNA pipeline R##
setwd("~/Dropbox/Horse_Transcriptome/downloads")
overallExpression <- read.table("dataSummary", header=T, stringsAsFactors=F)
overallExpression$transcriptName=rownames(overallExpression)
tissue_specific_exp <- read.table("backmapping_stats/allTissues_isoformTPM", header=T, stringsAsFactors=F)
#remove mt entries
tissue_specific_exp <- tissue_specific_exp[-c(1,2),]
rownames(tissue_specific_exp) <- c()
#load the different transcript categories
novel_I <- read.table("novelAnn/sup/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.sup", header=T, stringsAsFactors=F)
novel_II <- read.table("novelAnn/unsup.cons/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.unsup.cons", header=T, stringsAsFactors=F)
novel_III <- read.table("novelAnn/unsup.uncons.ORF/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.unsup.uncons.ORF", header=T, stringsAsFactors=F)
intergenic <- read.table("unsup", header=F, stringsAsFactors=F)

#attach global gene expression values and lenth 
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

#filter out transcripts GLOBALLY by <200 nt and with TPM < 0.1
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
#must get exon numbers by merging with bed file

##V10 in bed is exon number
single_novel_I <- subset(novel_I_bed, c(length < 1000 & V10 <2))
single_novel_II <- subset(novel_II_bed, c(length < 1000 & V10 <2))
single_novel_III <- subset(novel_III_bed, c(length < 1000 & V10 <2))
single_intergenic <- subset(intergenic_bed, c(length < 1000 & V10 <2))
#Filter out expression of single-exon transcripts in a tissue-specific manner
single_novel_I <- single_novel_I[apply(single_novel_I,1,function(row) {any(row[c(4-11)] > 5)}),]
single_novel_II <- single_novel_II[apply(single_novel_II,1,function(row) {any(row[c(4-11)] > 5)}),]
single_novel_III <- single_novel_III[apply(single_novel_III,1,function(row) {any(row[c(4-11)] > 5)}),]
single_intergenic <- single_intergenic[apply(single_intergenic,1,function(row) {any(row[c(4-11)] > 5)}),]
#remove non-bed format columns and format properly
single_novel_I_bed <- single_novel_I[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(single_novel_I_bed)[4]<-paste("V4")
rownames(single_novel_I_bed) <- c()
single_novel_II_bed <- single_novel_II[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(single_novel_II_bed)[4]<-paste("V4")
rownames(single_novel_II_bed) <- c()
single_novel_III_bed <- single_novel_III[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(single_novel_III_bed)[4]<-paste("V4")
rownames(single_novel_III_bed) <- c()
single_intergenic_bed <- single_intergenic[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(single_intergenic_bed)[4]<-paste("V4")
rownames(single_intergenic_bed) <- c()

#reorder original bed files into proper bed format
novel_I_bed <- novel_I_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_I_bed)[4]<-paste("V4")
novel_II_bed <- novel_II_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_II_bed)[4]<-paste("V4")
novel_III_bed <- novel_III_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_III_bed)[4]<-paste("V4")
intergenic_bed <- intergenic_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(intergenic_bed)[4]<-paste("V4")

#remove this subset from bed file
require(dplyr)
novel_I_bed <- anti_join(novel_I_bed,single_novel_I_bed, by="V4")
novel_II_bed <- anti_join(novel_II_bed,single_novel_II_bed, by="V4")
novel_III_bed <- anti_join(novel_III_bed,single_novel_III_bed, by="V4")
intergenic_bed <- anti_join(intergenic_bed,single_intergenic_bed, by="V4")

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

