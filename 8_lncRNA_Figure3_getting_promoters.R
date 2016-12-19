setwd("~/lncRNA")
##obtaining lncRNA promoters < 1kb 5'
lncRNA_all_bed <- read.table("lncRNA_final.bed", header=F,stringsAsFactors=F )
lncRNA_all_bed ["up"] <- (lncRNA_all_bed [ ,2]-1000)
lncRNA_all_bed ["down"]<- (lncRNA_all_bed [ ,3] + 1000)
#convert all negative numbers into 0
lncRNA_all_bed$up[lncRNA_all_bed$up < 0] <- 0
#isolate the -1000 (5'UTR)
lncRNA_promoters.pos <- lncRNA_all_bed[(lncRNA_all_bed$V6=="+"),c("V1","up","V2","V4","V5","V6","V7","V8","V9")]
lncRNA_promoters.neg <- lncRNA_all_bed[(lncRNA_all_bed$V6=="-"),c("V1","V3","down","V4","V5","V6","V7","V8","V9")]
names(lncRNA_promoters.neg)=names(lncRNA_promoters.pos)
lncRNA_promoters=rbind(lncRNA_promoters.pos,lncRNA_promoters.neg)
lncRNA_promoters <- lncRNA_promoters[with(lncRNA_promoters, order(V1, up)), ]

#isolate the +1000 (3'UTR)
lncRNA_down.pos <- lncRNA_all_bed[(lncRNA_all_bed$V6=="+"),c("V1","V3","down","V4","V5","V6","V7","V8","V9")]
lncRNA_down.neg <- lncRNA_all_bed[(lncRNA_all_bed$V6=="-"),c("V1","up","V2","V4","V5","V6","V7","V8","V9")]
names(lncRNA_down.neg)=names(lncRNA_down.pos)
lncRNA_down=rbind(lncRNA_down.pos,lncRNA_down.neg)
lncRNA_down <- lncRNA_down[with(lncRNA_down, order(V1, V3)), ]

write.table(lncRNA_promoters, "lncRNA_promoters.bed", row.names=F, col.names=F, quote=F,, sep = "\t")
write.table(lncRNA_down, "lncRNA_down.bed", row.names=F, col.names=F, quote=F,, sep = "\t")

#getting extended version of lncRNA bed in both directions
lncRNA_bed_extended <- lncRNA_all_bed[ ,c("V1","up","down","V4","V5","V6","V7","V8","V9","V10","V11","V12")]
names(lncRNA_bed_extended)[2]<-paste("V2")
names(lncRNA_bed_extended)[3]<-paste("V3")
#convert all negative numbers into 0
lncRNA_bed_extended$V2[lncRNA_bed_extended$V2 < 0] <- 0
#order chr properly
lncRNA_bed_extended <- lncRNA_bed_extended[with(lncRNA_bed_extended, order(V1, V2)), ]
write.table(lncRNA_bed_extended, "lncRNA_bed_extended.bed", row.names=F, col.names=F, quote=F, sep = "\t")


## create bed files of coding transcripts
candidateNonCoding_ids <- read.table("candidateNonCoding_ids", header=F,stringsAsFactors=F )
refined_bed <- read.table("inputs/allTissues_BED/refined.bed", header=F, stringsAsFactors=F)
refined_codingRNA_bed <- refined_bed[!(refined_bed$V4 %in% candidateNonCoding_ids$V1),]
refined_codingRNA_bed <- refined_codingRNA_bed[with(refined_codingRNA_bed, order(V1, V2)), ]
write.table(refined_codingRNA_bed,"refined_codingRNA.bed",row.names=F, col.names=F, quote=F, sep = "\t")

refined_codingRNA_bed["minus"] <- (refined_codingRNA_bed[ ,2]-1000)
refined_codingRNA_bed["add"]<- (refined_codingRNA_bed[ ,3] + 1000)
#Now to isolate the -1000 (5'UTR)
refined_codingRNA_bed_5.pos <- refined_codingRNA_bed[(refined_codingRNA_bed$V6=="+"),c("V1","minus","V2","V4","V5","V6","V7","V8","V9")]
refined_codingRNA_bed_5.neg <- refined_codingRNA_bed[(refined_codingRNA_bed$V6=="-"),c("V1","V3","add","V4","V5","V6","V7","V8","V9")]
names(refined_codingRNA_bed_5.neg)=names(refined_codingRNA_bed_5.pos)
refined_codingRNA_bed_5=rbind(refined_codingRNA_bed_5.pos,refined_codingRNA_bed_5.neg)
refined_codingRNA_bed_5$minus[refined_codingRNA_bed_5$minus < 0] <- 0
#order chr properly
refined_codingRNA_bed_5 <- refined_codingRNA_bed_5[with(refined_codingRNA_bed_5, order(V1, minus)), ]
write.table(refined_codingRNA_bed_5, "refined_codingRNA_5.bed", row.names=F, col.names=F, quote=F, sep = "\t")
#Now to isolate the +1000 (3'UTR)
refined_codingRNA_bed_3.pos <- refined_codingRNA_bed[(refined_codingRNA_bed$V6=="+"),c("V1","V3","add","V4","V5","V6","V7","V8","V9")]
refined_codingRNA_bed_3.neg <- refined_codingRNA_bed[(refined_codingRNA_bed$V6=="-"),c("V1","minus","V2","V4","V5","V6","V7","V8","V9")]
names(refined_codingRNA_bed_3.neg)=names(refined_codingRNA_bed_3.pos)
refined_codingRNA_bed_3=rbind(refined_codingRNA_bed_3.pos,refined_codingRNA_bed_3.neg)
refined_codingRNA_bed_3$V3[refined_codingRNA_bed_3$V3 < 0] <- 0
#order chr properly
refined_codingRNA_bed_3 <- refined_codingRNA_bed_3[with(refined_codingRNA_bed_3, order(V1, V3)), ]
write.table(refined_codingRNA_bed_3, "refined_codingRNA_3.bed", row.names=F, col.names=F, quote=F, sep = "\t")

