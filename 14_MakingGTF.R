##making GTF for the whole transcriptome
setwd("~/Dropbox/lncRNA")
lncRNA_gtf <- read.table("lncRNA_new.gtf", header=F, stringsAsFactors=F)
refined_nolncRNA_gtf  <- read.table("refined_nolncRNA.gtf", header=F, stringsAsFactors=F)
lncRNA_gtf["V2"] <- "lncRNA"
refined_nolncRNA_gtf ["V2"] <- "PCT"

total_transcriptome <- rbind(lncRNA_gtf,refined_nolncRNA_gtf)
#order chr properly
total_transcriptome[, c("V4")] <- sapply(total_transcriptome[, c("V4")], as.numeric)
total_transcriptome <- total_transcriptome[with(total_transcriptome, order(V1, V4)), ]
write.table(total_transcriptome, "ecTranscriptome.gtf", row.names=F, col.names=F, quote=F, sep = "\t")
