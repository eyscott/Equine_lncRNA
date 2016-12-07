setwd("~/lncRNA")
require(tidyr)
require(dplyr)
require(stringr)
P_blastn <- read.table("lncRNA_hg_nt.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
P_blastn <- P_blastn[!duplicated(P_blastn),]
names(P_blastn)=c("V5","sseqid","pident","qcovs","evalue")
P_noDups_id <- read.table("P_id.bed", header=F, stringsAsFactors=F,sep = "\t")
P_blastn_All <- merge(P_blastn,P_noDups_id, by="V5")

P_blastn_bed <- P_blastn_All[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13")]

P_lncRNA <- subset(P_blastn_All, qcovs >20)#makes a threshold of 20% coverage
P_lncRNA <- P_lncRNA[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13")] #amounts to 264 transcripts

P_final_bed <- anti_join(P_noDups_id,P_lncRNA,by=c("V1","V2","V3","V4"))

write.table(P_blastn_bed, "P_lncRNAall.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(P_lncRNA, "P_lncRNA.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(P_final_bed, "P_final.bed", row.names=F, col.names=F, quote=F, sep = "\t")
