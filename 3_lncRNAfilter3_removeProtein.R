require(tidyr)
require(dplyr)
require(stringr)
###Looking at BLASTp results
##looking at blastp.outfmt6 tables
novel_I_blastp <- read.table("novel_I_sprot.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
novel_II_blastp <- read.table("novel_II_sprot.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
novel_III_blastp <- read.table("novel_III_sprot.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
intergenic_blastp <- read.table("intergenic_sprot.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
known_lncRNA_blastp <- read.table("known_ncRNA_sprot.outfmt6", header=F, stringsAsFactors=F,sep = "\t")

#parse out column V1 to get gene IDs
novel_I_blastp_temp1=data.frame(str_split_fixed(novel_I_blastp$V1, ":", 4))
novel_I_blastp_bed <- cbind(novel_I_blastp_temp1$X3,novel_I_blastp[,c("V2","V11")])
names(novel_I_blastp_bed)=c("V4","target","sig")

novel_II_blastp_temp1=data.frame(str_split_fixed(novel_II_blastp$V1, ":", 4))
novel_II_blastp_bed <- cbind(novel_II_blastp_temp1$X3,novel_II_blastp[,c("V2","V11")])
names(novel_II_blastp_bed)=c("V4","target","sig")

novel_III_blastp_temp1=data.frame(str_split_fixed(novel_III_blastp$V1, ":", 4))
novel_III_blastp_bed <- cbind(novel_III_blastp_temp1$X3,novel_III_blastp[,c("V2","V11")])
names(novel_III_blastp_bed)=c("V4","target","sig")

intergenic_blastp_temp1=data.frame(str_split_fixed(intergenic_blastp$V1, ":", 4))
intergenic_blastp_bed <- cbind(intergenic_blastp_temp1$X3,intergenic_blastp[,c("V2","V11")])
names(intergenic_blastp_bed)=c("V4","target","sig")

known_lncRNA_blastp_temp1=data.frame(str_split_fixed(known_lncRNA_blastp$V1, ":", 4))
known_lncRNA_blastp_bed <- cbind(known_lncRNA_blastp_temp1$X3,known_lncRNA_blastp[,c("V2","V11")])
names(known_lncRNA_blastp_bed)=c("V4","target","sig")

#write blastp tables into comparable format with hmmersearch results
write.table(novel_I_blastp_bed, "novel_I_blastp.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_blastp_bed, "novel_II_blastp.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_blastp_bed, "novel_III_blastp.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_blastp_bed, "intergenic_blastp.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(known_lncRNA_blastp_bed, "known_ncRNA_blastp.bed", row.names=F, col.names=F, quote=F, sep = "\t")

##look at HMMER PFAM results
novel_I_pfam <- read.table("novel_I_pfam_new.tblout", header=F, stringsAsFactors=F)
novel_II_pfam <- read.table("novel_II_pfam_new.tblout", header=F, stringsAsFactors=F)
novel_III_pfam <- read.table("novel_III_pfam_new.tblout", header=F, stringsAsFactors=F)
intergenic_pfam <- read.table("intergenic_pfam_new.tblout", header=F, stringsAsFactors=F)
known_lncRNA_pfam <- read.table("known_ncRNA_pfam_new.tblout", header=F, stringsAsFactors=F)

###processing hmmsearch .tblout output to be compatible with blastp
novel_I_pfam_temp1=data.frame(str_split_fixed(novel_I_pfam$V1, ":", 4))
novel_I_pfam_bed <- cbind(novel_I_pfam_temp1$X3,novel_I_pfam[,c("V3","V5")])
names(novel_I_pfam_bed)=c("V4","target","sig")

novel_II_pfam_temp1=data.frame(str_split_fixed(novel_II_pfam$V1, ":", 4))
novel_II_pfam_bed <- cbind(novel_II_pfam_temp1$X3,novel_II_pfam[,c("V3","V5")])
names(novel_II_pfam_bed)=c("V4","target","sig")

novel_III_pfam_temp1=data.frame(str_split_fixed(novel_III_pfam$V1, ":", 4))
novel_III_pfam_bed <- cbind(novel_III_pfam_temp1$X3,novel_III_pfam[,c("V3","V5")])
names(novel_III_pfam_bed)=c("V4","target","sig")

intergenic_pfam_temp1=data.frame(str_split_fixed(intergenic_pfam$V1, ":", 4))
intergenic_pfam_bed <- cbind(intergenic_pfam_temp1$X3,intergenic_pfam[,c("V3","V5")])
names(intergenic_pfam_bed)=c("V4","target","sig")

known_lncRNA_pfam_temp1=data.frame(str_split_fixed(known_lncRNA_pfam$V1, ":", 4))
known_lncRNA_pfam_bed <- cbind(known_lncRNA_pfam_temp1$X3,known_lncRNA_pfam[,c("V3","V5")])
names(known_lncRNA_pfam_bed)=c("V4","target","sig")

write.table(novel_I_pfam_bed, "novel_I_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_pfam_bed, "novel_II_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_pfam_bed, "novel_III_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_pfam_bed, "intergenic_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(known_lncRNA_pfam_bed, "known_ncRNA_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")

##looking at blastn.outfmt6 tables
novel_I_blastn <- read.table("novel_I_hg38_cdna.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
novel_II_blastn <- read.table("novel_II_hg38_cdna.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
novel_III_blastn <- read.table("novel_III_hg38_cdna.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
intergenic_blastn <- read.table("intergenic_hg38_cdna.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
known_lncRNA_blastn <- read.table("known_ncRNA_hg38_cdna.outfmt6", header=F, stringsAsFactors=F,sep = "\t")

#parse out column V1 to get gene IDs
novel_I_blastn_temp1=data.frame(str_split_fixed(novel_I_blastn$V1, ":", 4))
novel_I_blastn_bed <- cbind(novel_I_blastn_temp1$X3,novel_I_blastn[,c("V2","V11")])
names(novel_I_blastn_bed)=c("V4","target","sig")

novel_II_blastn_temp1=data.frame(str_split_fixed(novel_II_blastn$V1, ":", 4))
novel_II_blastn_bed <- cbind(novel_II_blastn_temp1$X3,novel_II_blastn[,c("V2","V11")])
names(novel_II_blastn_bed)=c("V4","target","sig")

novel_III_blastn_temp1=data.frame(str_split_fixed(novel_III_blastn$V1, ":", 4))
novel_III_blastn_bed <- cbind(novel_III_blastn_temp1$X3,novel_III_blastn[,c("V2","V11")])
names(novel_III_blastn_bed)=c("V4","target","sig")

intergenic_blastn_temp1=data.frame(str_split_fixed(intergenic_blastn$V1, ":", 4))
intergenic_blastn_bed <- cbind(intergenic_blastn_temp1$X3,intergenic_blastn[,c("V2","V11")])
names(intergenic_blastn_bed)=c("V4","target","sig")

known_lncRNA_blastn_temp1=data.frame(str_split_fixed(known_lncRNA_blastn$V1, ":", 4))
known_lncRNA_blastn_bed <- cbind(known_lncRNA_blastn_temp1$X3,known_lncRNA_blastn[,c("V2","V11")])
names(known_lncRNA_blastn_bed)=c("V4","target","sig")

write.table(novel_I_blastn_bed, "novel_I_blastn.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_blastn_bed, "novel_II_blastn.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_blastn_bed, "novel_III_blastn.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_blastn_bed, "intergenic_blastn.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(known_lncRNA_blastn_bed, "known_ncRNA_blastn.bed", row.names=F, col.names=F, quote=F, sep = "\t")

#read in the remaining transcripts after filter 2
novel_I_bed_f2=read.table("novel_I_f2.bed", header=F, colClasses = "character",sep = "\t")
novel_II_bed_f2=read.table("novel_II_f2.bed", header=F, colClasses = "character",sep = "\t")
novel_III_bed_f2=read.table("novel_III_f2.bed", header=F, colClasses = "character",sep = "\t")
intergenic_bed_f2=read.table("intergenic_f2.bed", header=F, colClasses = "character",sep = "\t")
known_lncRNA_bed_f2=read.table("known_ncRNA_f2.bed", header=F, colClasses = "character",sep = "\t")

#combine all
all_f2_lncRNA <- rbind(novel_I_bed_f2,novel_II_bed_f2,novel_III_bed_f2,intergenic_bed_f2,known_lncRNA_bed_f2)
all_f2_lncRNA = all_f2_lncRNA[!duplicated(all_f2_lncRNA),]

##merge the protein findings
novel_I_P <- rbind(novel_I_blastp_bed,novel_I_pfam_bed,novel_I_blastn_bed)
novel_II_P <- rbind(novel_II_blastp_bed,novel_II_pfam_bed,novel_II_blastn_bed)
novel_III_P <- rbind(novel_III_blastp_bed,novel_III_pfam_bed,novel_III_blastn_bed)
intergenic_P <- rbind(intergenic_blastp_bed,intergenic_pfam_bed,intergenic_blastn_bed)
known_lncRNA_P <- rbind(known_lncRNA_blastp_bed,known_lncRNA_pfam_bed,known_lncRNA_blastn_bed)

######################
#removing duplicates from the two protein searches because they are causing
#problems with anti_join, as is the order of the chrs
novel_I_P_noDups <- novel_I_P[!duplicated(novel_I_P$V4),]
novel_II_P_noDups <- novel_II_P[!duplicated(novel_II_P$V4),]
novel_III_P_noDups <- novel_III_P[!duplicated(novel_III_P$V4),]
intergenic_P_noDups <- intergenic_P[!duplicated(intergenic_P$V4),]
known_lncRNA_P_noDups <- known_lncRNA_P[!duplicated(known_lncRNA_P$V4),]
###########
#merge with TCONS names
novel_I_P_noDups_bed <- merge(novel_I_P_noDups,all_f2_lncRNA,by="V4")
novel_I_P_noDups_bed = novel_I_P_noDups_bed[,colnames(all_f2_lncRNA)]
novel_II_P_noDups_bed <- merge(novel_II_P_noDups,all_f2_lncRNA,by="V4")
novel_II_P_noDups_bed = novel_II_P_noDups_bed[,colnames(all_f2_lncRNA)]
novel_III_P_noDups_bed <- merge(novel_III_P_noDups,all_f2_lncRNA,by="V4")
novel_III_P_noDups_bed = novel_III_P_noDups_bed[,colnames(all_f2_lncRNA)]
intergenic_P_noDups_bed <- merge(intergenic_P_noDups,all_f2_lncRNA,by="V4")
intergenic_P_noDups_bed = intergenic_P_noDups_bed[,colnames(all_f2_lncRNA)]
known_lncRNA_P_noDups_bed <- merge(known_lncRNA_P_noDups,all_f2_lncRNA,by="V4")
known_lncRNA_P_noDups_bed = known_lncRNA_P_noDups_bed[,colnames(all_f2_lncRNA)]

P_noDups_id <-rbind(data.frame(id="novel_I",novel_I_P_noDups_bed),
                 data.frame(id="novel_II",novel_II_P_noDups_bed),
                 data.frame(id="novel_III",novel_III_P_noDups_bed),
                 data.frame(id="intergenic",intergenic_P_noDups_bed),
                 data.frame(id="known",known_lncRNA_P_noDups_bed))

P_noDups <- data.frame(P_noDups_id[ ,-1])
P_noDups = P_noDups[!duplicated(P_noDups),]
P_noDups[, c("V2")] <- sapply(P_noDups[, c("V2")], as.numeric)
P_noDups <- P_noDups[with(P_noDups, order(V1,V2)), ]

write.table(novel_I_P_noDups, "novel_I_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_P_noDups, "novel_II_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_P_noDups, "novel_III_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_P_noDups, "intergenic_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(known_lncRNA_P_noDups, "known_ncRNA_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(P_noDups, "P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(P_noDups_id, "P_id.bed", row.names=F, col.names=F, quote=F, sep = "\t")


#performing anti_join to get rid of any transcript which had a hit in blastp or hmmsearch && make the f4 outputs into bed files
novel_I_bed <- anti_join(novel_I_bed_f2,novel_I_P_noDups, by="V4")
novel_II_bed <- anti_join(novel_II_bed_f2,novel_II_P_noDups, by="V4")
novel_III_bed <- anti_join(novel_III_bed_f2,novel_III_P_noDups, by="V4")
intergenic_bed <- anti_join(intergenic_bed_f2,intergenic_P_noDups, by="V4")
known_lncRNA_bed <- anti_join(known_lncRNA_bed_f2,known_lncRNA_P_noDups, by="V4")

novel_I_bed[, c("V2")] <- sapply(novel_I_bed[, c("V2")], as.numeric)
novel_II_bed[, c("V2")] <- sapply(novel_II_bed[, c("V2")], as.numeric)
novel_III_bed[, c("V2")] <- sapply(novel_III_bed[, c("V2")], as.numeric)
intergenic_bed[, c("V2")] <- sapply(intergenic_bed[, c("V2")], as.numeric)
known_lncRNA_bed[, c("V2")] <- sapply(known_lncRNA_bed[, c("V2")], as.numeric)

novel_I_bed <- novel_I_bed[with(novel_I_bed, order(V1, V2)), ]
novel_II_bed <- novel_II_bed[with(novel_II_bed, order(V1, V2)), ]
novel_III_bed <- novel_III_bed[with(novel_III_bed, order(V1, V2)), ]
intergenic_bed <- intergenic_bed[with(intergenic_bed, order(V1, V2)), ]
known_lncRNA_bed <- known_lncRNA_bed[with(known_lncRNA_bed, order(V1, V2)), ]

write.table(novel_I_bed, "novel_I_f3.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_bed, "novel_II_f3.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_bed, "novel_III_f3.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_bed, "intergenic_f3.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(known_lncRNA_bed, "known_ncRNA_f3.bed", row.names=F, col.names=F, quote=F, sep = "\t")


#now lets concat. them with labels to make them easier to make figures with later
lncRNA_all_Cat <-rbind(data.frame(id="novel_I",novel_I_bed),
                       data.frame(id="novel_II",novel_II_bed),
                       data.frame(id="novel_III",novel_III_bed),
                       data.frame(id="intergenic",intergenic_bed),
                       data.frame(id="known",known_lncRNA_bed))
                       
lncRNA_all_Cat <- lncRNA_all_Cat[with(lncRNA_all_Cat, order(V2, V3)), ]
write.table(lncRNA_all_Cat, "lncRNA_f3_IDs", row.names=F, col.names=F, quote=F, sep = "\t")
##now just making a table of f3, sub-divided into the protein coding vs non protein coding found in filter 3
all_ID <-rbind(data.frame(id="novel_I_lncRNA",novel_I_bed),
               data.frame(id="novel_I_genes",novel_I_P_noDups_bed),
               data.frame(id="novel_II_lncRNA",novel_II_bed),
               data.frame(id="novel_II_genes",novel_II_P_noDups_bed),
               data.frame(id="novel_III_lncRNA",novel_III_bed),
               data.frame(id="novel_III_genes",novel_III_P_noDups_bed),
               data.frame(id="intergenic_lncRNA",intergenic_bed),
               data.frame(id="intergenic_genes",intergenic_P_noDups_bed),
               data.frame(id="known_lncRNA",known_lncRNA_bed),
               data.frame(id="known_genes",known_lncRNA_P_noDups_bed))
               
write.table(all_ID, "all_cats_PandnoP", row.names=F, col.names=F, quote=F, sep = "\t")
