##lncRNA pipeline R##
setwd("~/Desktop/lncRNA")
overallExpression <- read.table("dataSummary", header=T, stringsAsFactors=F)

#load the different transcript categories
novel_I <- read.table("novel_I_genes.txt", header=F, stringsAsFactors=F)
novel_II <- read.table("novel_II_genes.txt", header=F, stringsAsFactors=F)
novel_III <- read.table("novel_III_genes.txt", header=F, stringsAsFactors=F)
intergenic <- 
  
#attach expression values and gene lengths
novel_I <- merge(overallExpression, novel_I, by.x="transcriptName",by.y="V2" )
novel_II <- merge(overallExpression, novel_II, by.x="transcriptName",by.y="V2" )
novel_III <- merge(overallExpression, novel_III, by.x="transcriptName",by.y="V2" )
intergenic

#slim down to columns you want
keeps <- c("transcriptName","length", "calcTPM")
novel_I <- unique(novel_I[keeps])
novel_II <- unique(novel_II[keeps])
novel_III <- unique(novel_III[keeps])
intergenic

#filter out transcripts <200 nt and with TPM < 0.1
novel_I_f1 <- subset(novel_I, c(calcTPM > 0.1 & length > 200))
novel_II_f1 <- subset(novel_II, c(calcTPM > 0.1 & length > 200))
novel_III_f1 <- subset(novel_III, c(calcTPM > 0.1 & length > 200))
intergenic_f1
#Make each into a bed file to use with bedtools
unfiltered_bed <- read.table("unfiltered_Alltissues_Assembly.bed", header=F, stringsAsFactors=F)
novel_I_bed <- merge(novel_I_f1, unfiltered_bed, by.x="transcriptName",by.y="V4" )
novel_II_bed <- merge(novel_II_f1, unfiltered_bed, by.x="transcriptName",by.y="V4" )
novel_III_bed <- merge(novel_III_f1, unfiltered_bed, by.x="transcriptName",by.y="V4" )

#reorder for proper bed format
novel_I_bed <- novel_I_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_I_bed)[4]<-paste("V4")
novel_II_bed <- novel_II_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_II_bed)[4]<-paste("V4")
novel_III_bed <- novel_III_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_III_bed)[4]<-paste("V4")
#order chr properly
novel_I_bed <- novel_I_bed[with(novel_I_bed, order(V1, V2)), ]
novel_II_bed <- novel_II_bed[with(novel_II_bed, order(V1, V2)), ]
novel_III_bed <- novel_III_bed[with(novel_III_bed, order(V1, V2)), ]

#write the bedfiles
write.table(novel_I_bed, "novel_I.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_bed, "novel_II.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_bed, "novel_III.bed", row.names=F, col.names=F, quote=F, sep = "\t")

#minus 1 kb and add 1 kb to end of each gene for merged transcriptome
merged_bed <- read.table("mergedTrans.bed", header=F,stringsAsFactors=F )
merged_bed["minus"] <- (merged_bed[ ,2]-1000)
merged_bed["add"]<- (merged_bed[ ,3] + 1000)
merged_bed_extended <- merged_bed[ ,c("V1","minus","add","V4","V5","V6","V7","V8","V9","V10","V11","V12")]
names(merged_bed_extended)[2]<-paste("V2")
names(merged_bed_extended)[3]<-paste("V3")
#convert all negative numbers into 0
merged_bed_extended$V2[merged_bed_extended$V2 < 0] <- 0
#order chr properly
merged_bed_extended <- merged_bed_extended[with(merged_bed_extended, order(V1, V2)), ]
write.table(merged_bed_extended, "merged_extended.bed", row.names=F, col.names=F, quote=F,, sep = "\t")
#this is to allow for removal of any transcripts within 1 kb of an annotated gene

#Now to isolate the -1000 (5'UTR)
merged_bed["minus"]<- (merged_bed[ ,2] - 1000)
merged_bed_5 <- merged_bed[ ,c("V1","minus","V2","V4","V5","V6","V7","V8","V9","V10","V11","V12")]
merged_bed_5$minus[merged_bed_5$minus < 0] <- 0
#order chr properly
merged_bed_5 <- merged_bed_5[with(merged_bed_5, order(V1, minus)), ]
write.table(merged_bed_5, "merged_5.bed", row.names=F, col.names=F, quote=F, sep = "\t")
#Now to isolate the +1000 (3'UTR)
merged_bed["add"]<- (merged_bed[ ,3] + 1000)
merged_bed_3 <- merged_bed[ ,c("V1","V3","add","V4","V5","V6","V7","V8","V9","V10","V11","V12")]
merged_bed_3 <- merged_bed_3[with(merged_bed_3, order(V1, V3)), ]
write.table(merged_bed_3, "merged_3.bed", row.names=F, col.names=F, quote=F, sep = "\t")
