setwd("~/lncRNA")
refined_bed <- read.table("inputs/allTissues_BED/refined.bed", header=F, stringsAsFactors=F)
f3 <- read.table("lncRNA_f3_IDs", header=F, stringsAsFactors=F)[-1]
#make sure the mergedtrans has none of the novel I,II,III lncRNA
library(dplyr)
names(f3) <- names(refined_bed)
refined_nolncRNA_bed <- anti_join(refined_bed,f3,by="V4")
#order chr properly
refined_nolncRNA_bed <- refined_nolncRNA_bed[with(refined_nolncRNA_bed, order(V1, V2)), ]
write.table(refined_nolncRNA_bed,"refined_nolncRNA.bed",row.names=F, col.names=F, quote=F, sep = "\t")

##observing overlap of lncRNA with transcriptome...done in R
#minus 1 kb and add 1 kb to end of each gene for merged transcriptome
refined_nolncRNA_bed["minus"] <- (refined_nolncRNA_bed[ ,2]-1000)
refined_nolncRNA_bed["add"]<- (refined_nolncRNA_bed[ ,3] + 1000)
refined_nolncRNA_bed_extended <- refined_nolncRNA_bed[ ,c("V1","minus","add","V4","V5","V6","V7","V8","V9","V10","V11","V12")]
names(refined_nolncRNA_bed_extended)[2]<-paste("V2")
names(refined_nolncRNA_bed_extended)[3]<-paste("V3")
#convert all negative numbers into 0
refined_nolncRNA_bed_extended$V2[refined_nolncRNA_bed_extended$V2 < 0] <- 0
#order chr properly
refined_nolncRNA_bed_extended <- refined_nolncRNA_bed_extended[with(refined_nolncRNA_bed_extended, order(V1, V2)), ]
write.table(refined_nolncRNA_bed_extended, "refined_nolncRNA_bed_extended.bed", row.names=F, col.names=F, quote=F, sep = "\t")
#this is to allow for removal of any transcripts within 1 kb of an annotated gene
#Now to isolate the -1000 (5'UTR)
refined_nolncRNA_bed_5 <- refined_nolncRNA_bed[ ,c("V1","minus","V2","V4","V5","V6","V7","V8","V9","V10","V11","V12")]
refined_nolncRNA_bed_5$minus[refined_nolncRNA_bed_5$minus < 0] <- 0
#order chr properly
refined_nolncRNA_bed_5 <- refined_nolncRNA_bed_5[with(refined_nolncRNA_bed_5, order(V1, minus)), ]
write.table(refined_nolncRNA_bed_5, "refined_nolncRNA_5.bed", row.names=F, col.names=F, quote=F, sep = "\t")
#Now to isolate the +1000 (3'UTR)
refined_nolncRNA_bed_3 <- refined_nolncRNA_bed[ ,c("V1","V3","add","V4","V5","V6","V7","V8","V9","V10","V11","V12")]
refined_nolncRNA_bed_3 <- refined_nolncRNA_bed_3[with(refined_nolncRNA_bed_3, order(V1, V3)), ]
write.table(refined_nolncRNA_bed_3, "refined_nolncRNA_3.bed", row.names=F, col.names=F, quote=F, sep = "\t")
