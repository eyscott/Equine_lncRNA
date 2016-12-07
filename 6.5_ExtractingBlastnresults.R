setwd("~/Desktop/lncRNA")
P_blastn <- read.table("lncRNA_hg_nt.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
#Get the protein that have hits to the lncRNA
P_blastn_temp1=data.frame(str_split_fixed(P_blastn$V1, ":", 2))
P_blastn_temp2=data.frame(str_split_fixed(P_blastn_temp1$X2, "\\(", 2))
P_blastn_temp3=data.frame(str_split_fixed(P_blastn_temp2$X1, "-", 2))
P_blastn_T_bed <- cbind(P_blastn_temp1$X1,P_blastn_temp3,P_blastn[,c("V2","V4")])
names(P_blastn_T_bed)=c("chr","start","stop","name","V4")
P_blastn_bed <- merge(P_blastn_T_bed,P_noDups_id, by=c("chr","start","stop"))
P_blastn_bed <- P_blastn_bed[ ,c("id","chr","start","stop","V4.y","V5","V6","V7","V8","V9","V10","V11","V12")]

P_bed <- read.table("P_id.bed", header=F, stringsAsFactors=F,sep = "\t")
names(P_bed) <- c("chr","start","stop","id","V4","V5","V6","V7","V8","V9","V10","V11","V12")

P_lncRNA <- subset(P_blastn_T_bed, V4 >20)#makes a threshold of 20% coverage
P_lncRNA <- merge(P_lncRNA,P_bed, by=c("chr","start","stop"))
P_lncRNA <- P_lncRNA[ ,c("id","chr","start","stop","V4.y","V5","V6","V7","V8","V9","V10","V11","V12")]
P_lncRNA <- P_lncRNA[!duplicated(P_lncRNA$V4.y),] #amounts to 264 transcripts
write.table(P_blastn_bed, "P_lncRNAall.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(P_lncRNA, "P_lncRNA.bed", row.names=F, col.names=F, quote=F, sep = "\t")
