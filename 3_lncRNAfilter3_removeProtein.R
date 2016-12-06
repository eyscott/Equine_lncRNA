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
novel_I_blastp_new <- separate(data = novel_I_blastp, col = V1, into = c("gene","chrs","g","m"), sep = "::")
novel_I_blastp_new <- novel_I_blastp_new[ ,c("chrs","V4","V2")] 
novel_I_blastp_new <- novel_I_blastp_new[!duplicated(novel_I_blastp_new),]
novel_I_V1 <- data.frame(gsub('.{3}$', '', novel_I_blastp_new$chrs))
novel_I_V1new <- separate(novel_I_V1, gsub....3.........novel_I_blastp_new.chrs.,into=c("start","stop"), sep="-")
novel_I_V1new <- separate(novel_I_V1new, start,into=c("chr","start"), sep=":")
novel_I_V2 <- separate(novel_I_blastp_new,col = V2, into=c("sp","Q","name"), sep=c(2,10))
novel_I_blastp_bed <- cbind(novel_I_blastp_new,novel_I_V1new,novel_I_V2)
novel_I_blastp_bed <- novel_I_blastp_bed[ ,c("chr","start","stop","name","V4")] 

novel_II_blastp_new <- separate(data = novel_II_blastp, col = V1, into = c("gene","chrs","g","m"), sep = "::")
novel_II_blastp_new <- novel_II_blastp_new[ ,c("chrs","V4","V2")] 
novel_II_blastp_new <- novel_II_blastp_new[!duplicated(novel_II_blastp_new),]
novel_II_V1 <- data.frame(gsub('.{3}$', '', novel_II_blastp_new$chr))
novel_II_V1new <- separate(novel_II_V1, gsub....3.........novel_II_blastp_new.chr.,into=c("start","stop"), sep="-")
novel_II_V1new <- separate(novel_II_V1new, start,into=c("chr","start"), sep=":")
novel_II_V2 <- separate(novel_II_blastp_new,col = V2, into=c("sp","Q","name"), sep=c(2,10))
novel_II_blastp_bed <- cbind(novel_II_blastp_new,novel_II_V1new,novel_II_V2)
novel_II_blastp_bed <- novel_II_blastp_bed[ ,c("chr","start","stop","name","V4")] 

novel_III_blastp_new <- separate(data = novel_III_blastp, col = V1, into = c("gene","chrs","g","m"), sep = "::")
novel_III_blastp_new <- novel_III_blastp_new[ ,c("chrs","V4","V2")] 
novel_III_blastp_new <- novel_III_blastp_new[!duplicated(novel_III_blastp_new),]
novel_III_V1 <- data.frame(gsub('.{3}$', '', novel_III_blastp_new$chr))
novel_III_V1new <- separate(novel_III_V1, gsub....3.........novel_III_blastp_new.chr.,into=c("start","stop"), sep="-")
novel_III_V1new <- separate(novel_III_V1new, start,into=c("chr","start"), sep=":")
novel_III_V2 <- separate(novel_III_blastp_new,col = V2, into=c("sp","Q","name"), sep=c(2,10))
novel_III_blastp_bed <- cbind(novel_III_blastp_new,novel_III_V1new,novel_III_V2)
novel_III_blastp_bed <- novel_III_blastp_bed[ ,c("chr","start","stop","name","V4")] 

intergenic_blastp_new <- separate(data = intergenic_blastp, col = V1, into = c("gene","chrs","g","m"), sep = "::")
intergenic_blastp_new <- intergenic_blastp_new[ ,c("chrs","V4","V2")]
intergenic_blastp_new <- intergenic_blastp_new[!duplicated(intergenic_blastp_new),]
intergenic_V1 <- data.frame(gsub('.{3}$', '', intergenic_blastp_new$chr))
intergenic_V1new <- separate(intergenic_V1, gsub....3.........intergenic_blastp_new.chr.,into=c("start","stop"), sep="-")
intergenic_V1new <- separate(intergenic_V1new, start,into=c("chr","start"), sep=":")
intergenic_V2 <- separate(intergenic_blastp_new,col = V2, into=c("sp","Q","name"), sep=c(2,10))
intergenic_blastp_bed <- cbind(intergenic_blastp_new,intergenic_V1new,intergenic_V2)
intergenic_blastp_bed <- intergenic_blastp_bed[ ,c("chr","start","stop","name","V4")] 

known_lncRNA_blastp_new <- separate(data = known_lncRNA_blastp, col = V1, into = c("gene","chrs","g","m"), sep = "::")
known_lncRNA_blastp_new <- known_lncRNA_blastp_new[ ,c("chrs","V4","V2")]
known_lncRNA_blastp_new <- known_lncRNA_blastp_new[!duplicated(known_lncRNA_blastp_new),]
known_lncRNA_V1 <- data.frame(gsub('.{3}$', '', known_lncRNA_blastp_new$chr))
known_lncRNA_V1new <- separate(known_lncRNA_V1, gsub....3.........known_lncRNA_blastp_new.chr.,into=c("start","stop"), sep="-")
known_lncRNA_V1new <- separate(known_lncRNA_V1new, start,into=c("chr","start"), sep=":")
known_lncRNA_V2 <- separate(known_lncRNA_blastp_new,col = V2, into=c("sp","Q","name"), sep=c(2,10))
known_lncRNA_blastp_bed <- cbind(known_lncRNA_blastp_new,known_lncRNA_V1new,known_lncRNA_V2)
known_lncRNA_blastp_bed <- known_lncRNA_blastp_bed[ ,c("chr","start","stop","name","V4")] 

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
novel_I_pfam_new <- separate(data = novel_I_pfam, col = V22, into = c("chr","starts","subj"), sep = ":")
novel_I_pfam_new <- novel_I_pfam_new[ ,c("chr","starts","V20","V5","V8", "V3")] 
novel_I_pfam_new <- novel_I_pfam_new[!duplicated(novel_I_pfam_new),]
novel_I_V22 <- data.frame(gsub('.{3}$', '', novel_I_pfam_new$start))
novel_I_V22new <- separate(novel_I_V22, gsub....3.........novel_I_pfam_new.start.,into=c("start","stop"), sep="-")
novel_I_length <- separate(data = novel_I_pfam_new, col = V20, into = c("len","length"), sep = ":")
novel_I_pfam_bed <- cbind(novel_I_pfam_new,novel_I_V22new, novel_I_length)
novel_I_pfam_bed <- novel_I_pfam_bed[ ,c("chr","start","stop","V3","length", "V5","V8")] 
novel_I_pfam_sub <- subset(novel_I_pfam_bed,c(V5 > 0.001 & V8 > 0.001))
novel_I_pfam_sub <- novel_I_pfam_sub[ ,c("chr","start","stop","V3","length")]

novel_II_pfam_new <- separate(data = novel_II_pfam, col = V22, into = c("chr","starts","subj"), sep = ":")
novel_II_pfam_new <- novel_II_pfam_new[ ,c("chr","starts","V20","V5","V8", "V3")] 
novel_II_pfam_new <- novel_II_pfam_new[!duplicated(novel_II_pfam_new),]
novel_II_V22 <- data.frame(gsub('.{3}$', '', novel_II_pfam_new$start))
novel_II_V22new <- separate(novel_II_V22, gsub....3.........novel_II_pfam_new.start.,into=c("start","stop"), sep="-")
novel_II_length <- separate(data = novel_II_pfam_new, col = V20, into = c("len","length"), sep = ":")
novel_II_pfam_bed <- cbind(novel_II_pfam_new,novel_II_V22new, novel_II_length)
novel_II_pfam_bed <- novel_II_pfam_bed[ ,c("chr","start","stop","V3","length", "V5","V8")] 
novel_II_pfam_sub <- subset(novel_II_pfam_bed,c(V5 > 0.001 & V8 > 0.001))
novel_II_pfam_sub <- novel_II_pfam_sub[ ,c("chr","start","stop","V3","length")]

novel_III_pfam_new <- separate(data = novel_III_pfam, col = V22, into = c("chr","starts","subj"), sep = ":")
novel_III_pfam_new <- novel_III_pfam_new[ ,c("chr","starts","V20","V5","V8", "V3")] 
novel_III_pfam_new <- novel_III_pfam_new[!duplicated(novel_III_pfam_new),]
novel_III_V22 <- data.frame(gsub('.{3}$', '', novel_III_pfam_new$start))
novel_III_V22new <- separate(novel_III_V22, gsub....3.........novel_III_pfam_new.start.,into=c("start","stop"), sep="-")
novel_III_length <- separate(data = novel_III_pfam_new, col = V20, into = c("len","length"), sep = ":")
novel_III_pfam_bed <- cbind(novel_III_pfam_new,novel_III_V22new, novel_III_length)
novel_III_pfam_bed <- novel_III_pfam_bed[ ,c("chr","start","stop","V3","length", "V5","V8")] 
novel_III_pfam_sub <- subset(novel_III_pfam_bed,c(V5 > 0.001 & V8 > 0.001))
novel_III_pfam_sub <- novel_III_pfam_sub[ ,c("chr","start","stop","V3","length")]

intergenic_pfam_new <- separate(data = intergenic_pfam, col = V22, into = c("chr","starts","subj"), sep = ":")
intergenic_pfam_new <- intergenic_pfam_new[ ,c("chr","starts","V20","V5","V8", "V3")] 
intergenic_pfam_new <- intergenic_pfam_new[!duplicated(intergenic_pfam_new),]
intergenic_V22 <- data.frame(gsub('.{3}$', '', intergenic_pfam_new$start))
intergenic_V22new <- separate(intergenic_V22, gsub....3.........intergenic_pfam_new.start.,into=c("start","stop"), sep="-")
intergenic_length <- separate(data = intergenic_pfam_new, col = V20, into = c("len","length"), sep = ":")
intergenic_pfam_bed <- cbind(intergenic_pfam_new,intergenic_V22new, intergenic_length)
intergenic_pfam_bed <- intergenic_pfam_bed[ ,c("chr","start","stop","V3","length", "V5","V8")] 
intergenic_pfam_sub <- subset(intergenic_pfam_bed,c(V5 > 0.001 & V8 > 0.001))
intergenic_pfam_sub <- intergenic_pfam_sub[ ,c("chr","start","stop","V3","length")]

known_lncRNA_pfam_new <- separate(data = known_lncRNA_pfam, col = V22, into = c("chr","starts","subj"), sep = ":")
known_lncRNA_pfam_new <- known_lncRNA_pfam_new[ ,c("chr","starts","V20","V5","V8", "V3")] 
known_lncRNA_pfam_new <- known_lncRNA_pfam_new[!duplicated(known_lncRNA_pfam),]
known_lncRNA_V22 <- data.frame(gsub('.{3}$', '', known_lncRNA_pfam_new$start))
known_lncRNA_V22new <- separate(known_lncRNA_V22, gsub....3.........known_lncRNA_pfam_new.start.,into=c("start","stop"), sep="-")
known_lncRNA_length <- separate(data = known_lncRNA_pfam_new, col = V20, into = c("len","length"), sep = ":")
known_lncRNA_pfam_bed <- cbind(known_lncRNA_pfam_new,known_lncRNA_V22new, known_lncRNA_length)
known_lncRNA_pfam_bed <- known_lncRNA_pfam_bed[ ,c("chr","start","stop","V3","length", "V5","V8")] 
known_lncRNA_pfam_sub <- subset(known_lncRNA_pfam_bed,c(V5 > 0.001 & V8 > 0.001))
known_lncRNA_pfam_sub <- known_lncRNA_pfam_sub[ ,c("chr","start","stop","V3","length")]


write.table(novel_I_pfam_sub, "novel_I_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_pfam_sub, "novel_II_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_pfam_sub, "novel_III_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_pfam_sub, "intergenic_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(known_lncRNA_pfam_sub, "known_ncRNA_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")

# order of table: chr,starts,stop,V3-Name,length

##looking at blastn.outfmt6 tables
novel_I_blastn <- read.table("novel_I_hg38_cdna.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
novel_II_blastn <- read.table("novel_II_hg38_cdna.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
novel_III_blastn <- read.table("novel_III_hg38_cdna.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
intergenic_blastn <- read.table("intergenic_hg38_cdna.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
known_lncRNA_blastn <- read.table("known_ncRNA_hg38_cdna.outfmt6", header=F, stringsAsFactors=F,sep = "\t")

#parse out column V1 to get gene IDs
novel_I_blastn_temp1=data.frame(str_split_fixed(novel_I_blastn$V1, ":", 2))
novel_I_blastn_temp2=data.frame(str_split_fixed(novel_I_blastn_temp1$X2, "\\(", 2))
novel_I_blastn_temp3=data.frame(str_split_fixed(novel_I_blastn_temp2$X1, "-", 2))
novel_I_blastn_bed <- cbind(novel_I_blastn_temp1$X1,novel_I_blastn_temp3,novel_I_blastn[,c("V2","V4")])
names(novel_I_blastn_bed)=c("chr","start","stop","name","V4")

novel_II_blastn_temp1=data.frame(str_split_fixed(novel_II_blastn$V1, ":", 2))
novel_II_blastn_temp2=data.frame(str_split_fixed(novel_II_blastn_temp1$X2, "\\(", 2))
novel_II_blastn_temp3=data.frame(str_split_fixed(novel_II_blastn_temp2$X1, "-", 2))
novel_II_blastn_bed <- cbind(novel_II_blastn_temp1$X1,novel_II_blastn_temp3,novel_II_blastn[,c("V2","V4")])
names(novel_II_blastn_bed)=c("chr","start","stop","name","V4")

novel_III_blastn_temp1=data.frame(str_split_fixed(novel_III_blastn$V1, ":", 2))
novel_III_blastn_temp2=data.frame(str_split_fixed(novel_III_blastn_temp1$X2, "\\(", 2))
novel_III_blastn_temp3=data.frame(str_split_fixed(novel_III_blastn_temp2$X1, "-", 2))
novel_III_blastn_bed <- cbind(novel_III_blastn_temp1$X1,novel_III_blastn_temp3,novel_III_blastn[,c("V2","V4")])
names(novel_III_blastn_bed)=c("chr","start","stop","name","V4")

intergenic_blastn_temp1=data.frame(str_split_fixed(intergenic_blastn$V1, ":", 2))
intergenic_blastn_temp2=data.frame(str_split_fixed(intergenic_blastn_temp1$X2, "\\(", 2))
intergenic_blastn_temp3=data.frame(str_split_fixed(intergenic_blastn_temp2$X1, "-", 2))
intergenic_blastn_bed <- cbind(intergenic_blastn_temp1$X1,intergenic_blastn_temp3,intergenic_blastn[,c("V2","V4")])
names(intergenic_blastn_bed)=c("chr","start","stop","name","V4")

known_lncRNA_blastn_temp1=data.frame(str_split_fixed(known_lncRNA_blastn$V1, ":", 2))
known_lncRNA_blastn_temp2=data.frame(str_split_fixed(known_lncRNA_blastn_temp1$X2, "\\(", 2))
known_lncRNA_blastn_temp3=data.frame(str_split_fixed(known_lncRNA_blastn_temp2$X1, "-", 2))
known_lncRNA_blastn_bed <- cbind(known_lncRNA_blastn_temp1$X1,known_lncRNA_blastn_temp3,known_lncRNA_blastn[,c("V2","V4")])
names(known_lncRNA_blastn_bed)=c("chr","start","stop","name","V4")


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
#all_f2_lncRNA <- rbind(novel_I_bed_f2,novel_II_bed_f2,novel_III_bed_f2,intergenic_bed_f2,known_lncRNA_bed_f2)
#format tables for comparison
trunc_headers=c("chr","start","stop")
names(novel_I_bed_f2)[1:3]=trunc_headers
names(novel_II_bed_f2)[1:3]=trunc_headers
names(novel_III_bed_f2)[1:3]=trunc_headers
names(intergenic_bed_f2)[1:3]=trunc_headers
names(known_lncRNA_bed_f2)[1:3]=trunc_headers
#names(all_f2_lncRNA)[1:3]=trunc_headers

novel_I_pfam_sub_trunc <-novel_I_pfam_sub[ ,trunc_headers]
novel_II_pfam_sub_trunc <-novel_II_pfam_sub[ ,trunc_headers]
novel_III_pfam_sub_trunc <-novel_III_pfam_sub[ ,trunc_headers]
intergenic_pfam_sub_trunc <-intergenic_pfam_sub[ ,trunc_headers]
known_lncRNA_pfam_sub_trunc <-known_lncRNA_pfam_sub[ ,trunc_headers]

novel_I_blastp_bed_trunc <-novel_I_blastp_bed[ ,trunc_headers]
novel_II_blastp_bed_trunc <-novel_II_blastp_bed[ ,trunc_headers]
novel_III_blastp_bed_trunc <-novel_III_blastp_bed[ ,trunc_headers]
intergenic_blastp_bed_trunc <-intergenic_blastp_bed[ ,trunc_headers]
known_lncRNA_blastp_bed_trunc <-known_lncRNA_blastp_bed[ ,trunc_headers]

novel_I_blastn_bed_trunc <-novel_I_blastn_bed[ ,trunc_headers]
novel_II_blastn_bed_trunc <-novel_II_blastn_bed[ ,trunc_headers]
novel_III_blastn_bed_trunc <-novel_III_blastn_bed[ ,trunc_headers]
intergenic_blastn_bed_trunc <- intergenic_blastn_bed[ ,trunc_headers]
known_lncRNA_blastn_bed_trunc <- known_lncRNA_blastn_bed[ ,trunc_headers]

##merge the protein findings
novel_I_P <- rbind(novel_I_blastp_bed_trunc,novel_I_pfam_sub_trunc,novel_I_blastn_bed_trunc)
novel_II_P <- rbind(novel_II_blastp_bed_trunc,novel_II_pfam_sub_trunc,novel_II_blastn_bed_trunc)
novel_III_P <- rbind(novel_III_blastp_bed_trunc,novel_III_pfam_sub_trunc,novel_III_blastn_bed_trunc)
intergenic_P <- rbind(intergenic_blastp_bed_trunc,intergenic_pfam_sub_trunc,intergenic_blastn_bed_trunc)
known_lncRNA_P <- rbind(known_lncRNA_blastp_bed_trunc,known_lncRNA_pfam_sub_trunc,known_lncRNA_blastn_bed_trunc)

######################
#removing duplicates from the two protein searches because they are causing
#problems with anti_join, as is the order of the chrs
novel_I_P_noDups <- novel_I_P[!duplicated(novel_I_P),]
novel_I_P_noDups <- novel_I_P_noDups[with(novel_I_P_noDups, order(chr,start)), ]

novel_II_P_noDups <- novel_II_P[!duplicated(novel_II_P),]
novel_II_P_noDups <- novel_II_P_noDups[with(novel_II_P_noDups, order(chr,start)), ]

novel_III_P_noDups <- novel_III_P[!duplicated(novel_III_P),]
novel_III_P_noDups <- novel_III_P_noDups[with(novel_III_P_noDups, order(chr,start)), ]

intergenic_P_noDups <- intergenic_P[!duplicated(intergenic_P),]
intergenic_P_noDups <- intergenic_P_noDups[with(intergenic_P_noDups, order(chr,start)), ]

known_lncRNA_P_noDups <- known_lncRNA_P[!duplicated(known_lncRNA_P),]
known_lncRNA_P_noDups <- known_lncRNA_P_noDups[with(known_lncRNA_P_noDups, order(chr,start)), ]

#P_noDups <-rbind(data.frame(id="novel_I",novel_I_P_noDups),
#                       data.frame(id="novel_II",novel_II_P_noDups),
#                       data.frame(id="novel_III",novel_III_P_noDups),
#                       data.frame(id="intergenic",intergenic_P_noDups),
#                       data.frame(id="known",known_lncRNA_P_noDups))
#merge with TCONS names
#P_noDups <- merge(P_noDups,all_f2_lncRNA,by=c("chr","start","stop"))
#P_noDups <- P_noDups[with(P_noDups, order(chr,start)), ]

write.table(novel_I_P_noDups, "novel_I_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_P_noDups, "novel_II_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_P_noDups, "novel_III_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_P_noDups, "intergenic_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(known_lncRNA_P_noDups, "known_ncRNA_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
#write.table(P_noDups, "P.bed", row.names=F, col.names=F, quote=F, sep = "\t")

#performing anti_join to get rid of any transcript which had a hit in blastp or hmmsearch && make the f4 outputs into bed files
novel_I_bed <- anti_join(novel_I_bed_f2,novel_I_P_noDups, by=c("chr","start","stop"))
novel_II_bed <- anti_join(novel_II_bed_f2,novel_II_P_noDups, by=c("chr","start","stop"))
novel_III_bed <- anti_join(novel_III_bed_f2,novel_III_P_noDups, by=c("chr","start","stop"))
intergenic_bed <- anti_join(intergenic_bed_f2,intergenic_P_noDups, by=c("chr","start","stop"))
known_lncRNA_bed <- anti_join(known_lncRNA_bed_f2,intergenic_P_noDups, by=c("chr","start","stop"))

novel_I_bed[, c("start")] <- sapply(novel_I_bed[, c("start")], as.numeric)
novel_II_bed[, c("start")] <- sapply(novel_II_bed[, c("start")], as.numeric)
novel_III_bed[, c("start")] <- sapply(novel_III_bed[, c("start")], as.numeric)
intergenic_bed[, c("start")] <- sapply(intergenic_bed[, c("start")], as.numeric)
known_lncRNA_bed[, c("start")] <- sapply(known_lncRNA_bed[, c("start")], as.numeric)

novel_I_bed <- novel_I_bed[with(novel_I_bed, order(chr, start)), ]
novel_II_bed <- novel_II_bed[with(novel_II_bed, order(chr, start)), ]
novel_III_bed <- novel_III_bed[with(novel_III_bed, order(chr, start)), ]
intergenic_bed <- intergenic_bed[with(intergenic_bed, order(chr, start)), ]
known_lncRNA_bed <- known_lncRNA_bed[with(known_lncRNA_bed, order(chr, start)), ]

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
                       
lncRNA_all_Cat <- lncRNA_all_Cat[with(lncRNA_all_Cat, order(chr, start)), ]
write.table(lncRNA_all_Cat, "lncRNA_f3_IDs", row.names=F, col.names=F, quote=F, sep = "\t")
##now just making a table of f3, sub-divided into the protein coding vs non protein coding found in filter 3
all_ID <-rbind(data.frame(id="novel_I_lncRNA",novel_I_bed[,1:3]),
               data.frame(id="novel_I_genes",novel_I_P_noDups),
               data.frame(id="novel_II_lncRNA",novel_II_bed[,1:3]),
               data.frame(id="novel_II_genes",novel_II_P_noDups),
               data.frame(id="novel_III_lncRNA",novel_III_bed[,1:3]),
               data.frame(id="novel_III_genes",novel_III_P_noDups),
               data.frame(id="intergenic_lncRNA",intergenic_bed[,1:3]),
               data.frame(id="intergenic_genes",intergenic_P_noDups),
               data.frame(id="known_lncRNA",known_lncRNA_bed[,1:3]),
               data.frame(id="known_genes",known_lncRNA_P_noDups))
               
write.table(all_ID, "all_cats_PandnoP", row.names=F, col.names=F, quote=F, sep = "\t")
