setwd("~/Dropbox/lncRNA")
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

#getting extended version of lncRNA bed in both directions
lncRNA_bed_extended <- lncRNA_all_bed[ ,c("V1","up","down","V4","V5","V6","V7","V8","V9","V10","V11","V12")]
names(lncRNA_bed_extended)[2]<-paste("V2")
names(lncRNA_bed_extended)[3]<-paste("V3")
#convert all negative numbers into 0
lncRNA_bed_extended$V2[lncRNA_bed_extended$V2 < 0] <- 0
#order chr properly
lncRNA_bed_extended <- lncRNA_bed_extended[with(lncRNA_bed_extended, order(V1, V2)), ]
setwd("~/Desktop/lncRNA")
write.table(lncRNA_bed_extended, "lncRNA_bed_extended.bed", row.names=F, col.names=F, quote=F, sep = "\t")

lncRNA_extend_bed <- read.table("lncRNA_bed_extended.bed", header=F) 
refined_extend_bed <- read.table("refined_nolncRNA_bed_extended.bed", header=F) 
equine_transcription.bed <- rbind(lncRNA_extend_bed,refined_extend_bed)
#order chr properly
equine_transcription.bed[, c("V2")] <- sapply(equine_transcription.bed[, c("V2")], as.numeric)
equine_transcription.bed <- equine_transcription.bed[with(equine_transcription.bed, order(V1, V2)), ]
write.table(equine_transcription.bed,"equine_transcription.bed",row.names=F, col.names=F, quote=F, sep = "\t")

ecu_chr <- read.table("equCab2.chrom.sizes.txt", header=T)
ecu_chr <- ecu_chr[with(ecu_chr, order(chr)), ]
write.table(ecu_chr,"equCab2.chrom.sizes",row.names=F, col.names=F, quote=F, sep = "\t")


