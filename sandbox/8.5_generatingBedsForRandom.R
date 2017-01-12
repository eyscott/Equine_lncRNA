setwd("~/lncRNA")
lncRNA_extend_bed <- read.table("lncRNA_bed_extended.bed", header=F) 
refined_extend_bed <- read.table("refined_nolncRNA_bed_extended.bed", header=F) 
equine_transcription.bed <- rbind(lncRNA_extend_bed,refined_extend_bed)
#order chr properly
equine_transcription.bed[, c("V2")] <- sapply(equine_transcription.bed[, c("V2")], as.numeric)
equine_transcription.bed <- equine_transcription.bed[with(equine_transcription.bed, order(V1, V2)), ]
write.table(equine_transcription.bed,"equine_transcription.bed",row.names=F, col.names=F, quote=F, sep = "\t")

ecu_chr <- read.table("inputs/equCab2.chrom.sizes", header=F, col.names=c("chr","size"))
ecu_chr <- ecu_chr[with(ecu_chr, order(chr)), ]
write.table(ecu_chr,"equCab2.chrom.sizes",row.names=F, col.names=F, quote=F, sep = "\t")


