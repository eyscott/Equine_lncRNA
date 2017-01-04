##making GTF from BED file
setwd("~/Dropbox/lncRNA")
lncRNA_bed <- read.table("lncRNA_final.bed", header=F, stringsAsFactors=F)
#make some new columns for GTF
lncRNA_bed["source"] <- "lncRNA"
lncRNA_bed["feature"] <- "transcript"
lncRNA_bed["score"] <- "."
lncRNA_bed["frame"] <- "."

# rearrange to make this into a GTF
lncRNA_gtf <- lncRNA_bed[ ,c("V1","source","feature","V2","V3","source","V6","frame","V4")]
write.table(lncRNA_gtf, "lncRNA.gtf", row.names=F, col.names=F, quote=F, sep = "\t")
