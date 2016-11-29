##observing overlap of lncRNA with transcriptome...done in R
#minus 1 kb and add 1 kb to end of each gene for merged transcriptome
setwd("~/Dropbox/Horse_Transcriptome/downloads")
merged_bed <- read.table("allTissues_BED/mergedTrans.bed", header=F,stringsAsFactors=F )
merged_bed["minus"] <- (merged_bed[ ,2]-1000)
merged_bed["add"]<- (merged_bed[ ,3] + 1000)
merged_bed_extended <- merged_bed[ ,c("V1","minus","add","V4","V5","V6","V7","V8","V9","V10","V11","V12")]
names(merged_bed_extended)[2]<-paste("V2")
names(merged_bed_extended)[3]<-paste("V3")
#convert all negative numbers into 0
merged_bed_extended$V2[merged_bed_extended$V2 < 0] <- 0
#order chr properly
merged_bed_extended <- merged_bed_extended[with(merged_bed_extended, order(V1, V2)), ]
write.table(merged_bed_extended, "merged_extended.bed", row.names=F, col.names=F, quote=F, sep = "\t")
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
