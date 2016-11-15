##lncRNA pipeline R##
setwd("~/Dropbox/Horse_Transcriptome/downloads")
overallExpression <- read.table("dataSummary", header=T, stringsAsFactors=F)
overallExpression$transcriptName=rownames(overallExpression)
tissue_specific_exp <- read.table("backmapping_stats/allTissues_isoformTPM", header=T, stringsAsFactors=F)
#remove mt entries
tissue_specific_exp <- tissue_specific_exp[-c(1,2),]
rownames(tissue_specific_exp) <- c()
#load the different transcript categories
novel_I <- read.table("novelAnn/sup/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.sup", header=T, stringsAsFactors=F)
novel_II <- read.table("novelAnn/unsup.cons/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.unsup.cons", header=T, stringsAsFactors=F)
novel_III <- read.table("novelAnn/unsup.uncons.ORF/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.unsup.uncons.ORF", header=T, stringsAsFactors=F)
intergenic <- read.table("unsup", header=F, stringsAsFactors=F)

#attach global gene expression values and lenth 
novel_I <- merge(overallExpression, novel_I, by.x="transcriptName",by.y="transcript.ID" )
novel_II <- merge(overallExpression, novel_II, by.x="transcriptName",by.y="transcript.ID" )
novel_III <- merge(overallExpression, novel_III, by.x="transcriptName",by.y="transcript.ID" )
intergenic <- merge(overallExpression, intergenic, by.x="transcriptName",by.y="V1" )

#slim down to columns you want
keeps <- c("transcriptName","length", "calcTPM")
novel_I <- unique(novel_I[keeps])
novel_II <- unique(novel_II[keeps])
novel_III <- unique(novel_III[keeps])
intergenic <- unique(intergenic[keeps])

#filter out transcripts GLOBALLY by <200 nt and with TPM < 0.1
novel_I_f1 <- subset(novel_I, c(calcTPM > 0.1 & length > 200))
novel_II_f1 <- subset(novel_II, c(calcTPM > 0.1 & length > 200))
novel_III_f1 <- subset(novel_III, c(calcTPM > 0.1 & length > 200))
intergenic_f1 <- subset(intergenic, c(calcTPM > 0.1 & length > 200))

#Make each into a bed file to use with bedtools
setwd("~/Dropbox/Horse_Transcriptome/downloads/allTissues_BED")
unfiltered_bed <- read.table("unfiltered_Alltissues_Assembly.bed", header=F, stringsAsFactors=F)
novel_I_bed <- merge(novel_I_f1, unfiltered_bed, by.x="transcriptName",by.y="V4" )
novel_II_bed <- merge(novel_II_f1, unfiltered_bed, by.x="transcriptName",by.y="V4" )
novel_III_bed <- merge(novel_III_f1, unfiltered_bed, by.x="transcriptName",by.y="V4" )
intergenic_bed <- merge(intergenic_f1, unfiltered_bed, by.x="transcriptName",by.y="V4" )

#applying more stringent criteria to single exon transcripts
#must get exon numbers by merging with bed file

##V10 in bed is exon number
single_novel_I <- subset(novel_I_bed, c(length < 1000 & V10 <2))
single_novel_II <- subset(novel_II_bed, c(length < 1000 & V10 <2))
single_novel_III <- subset(novel_III_bed, c(length < 1000 & V10 <2))
single_intergenic <- subset(intergenic_bed, c(length < 1000 & V10 <2))
#Filter out expression of single-exon transcripts in a tissue-specific manner
single_novel_I <- single_novel_I[apply(single_novel_I,1,function(row) {any(row[c(4-11)] > 5)}),]
single_novel_II <- single_novel_II[apply(single_novel_II,1,function(row) {any(row[c(4-11)] > 5)}),]
single_novel_III <- single_novel_III[apply(single_novel_III,1,function(row) {any(row[c(4-11)] > 5)}),]
single_intergenic <- single_intergenic[apply(single_intergenic,1,function(row) {any(row[c(4-11)] > 5)}),]
#remove non-bed format columns and format properly
single_novel_I_bed <- single_novel_I[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(single_novel_I_bed)[4]<-paste("V4")
rownames(single_novel_I_bed) <- c()
single_novel_II_bed <- single_novel_II[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(single_novel_II_bed)[4]<-paste("V4")
rownames(single_novel_II_bed) <- c()
single_novel_III_bed <- single_novel_III[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(single_novel_III_bed)[4]<-paste("V4")
rownames(single_novel_III_bed) <- c()
single_intergenic_bed <- single_intergenic[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(single_intergenic_bed)[4]<-paste("V4")
rownames(single_intergenic_bed) <- c()

#reorder original bed files into proper bed format
novel_I_bed <- novel_I_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_I_bed)[4]<-paste("V4")
novel_II_bed <- novel_II_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_II_bed)[4]<-paste("V4")
novel_III_bed <- novel_III_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_III_bed)[4]<-paste("V4")
intergenic_bed <- intergenic_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(intergenic_bed)[4]<-paste("V4")

#remove this subset from bed file
require(dplyr)
novel_I_bed <- anti_join(novel_I_bed,single_novel_I_bed, by="V4")
novel_II_bed <- anti_join(novel_II_bed,single_novel_II_bed, by="V4")
novel_III_bed <- anti_join(novel_III_bed,single_novel_III_bed, by="V4")
intergenic_bed <- anti_join(intergenic_bed,single_intergenic_bed, by="V4")

#order chr properly
novel_I_bed <- novel_I_bed[with(novel_I_bed, order(V1, V2)), ]
novel_II_bed <- novel_II_bed[with(novel_II_bed, order(V1, V2)), ]
novel_III_bed <- novel_III_bed[with(novel_III_bed, order(V1, V2)), ]
intergenic_bed <- intergenic_bed[with(intergenic_bed, order(V1, V2)), ]
unfiltered_bed <- unfiltered_bed[with(unfiltered_bed, order(V1, V2)), ]
#write the bedfiles
setwd("~/Desktop/lncRNA")
write.table(novel_I_bed, "novel_I.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_bed, "novel_II.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_bed, "novel_III.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_bed, "intergenic.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(unfiltered_bed, "unfiltered.bed", row.names=F, col.names=F, quote=F, sep = "\t")

##transcripts to 5' and 3'...reading the bedtools intersect files
setwd("~/Desktop/lncRNA")
novel_I_5 <- read.table("novel_I_5.bed", header=F, stringsAsFactors=F)
novel_I_3 <- read.table("novel_I_3.bed", header=F, stringsAsFactors=F)
novel_I_rest <- read.table("novel_I_ex.bed", header=F, stringsAsFactors=F)
novel_II_5 <- read.table("novel_II_5.bed", header=F, stringsAsFactors=F)
novel_II_3 <- read.table("novel_II_3.bed", header=F, stringsAsFactors=F)
novel_II_rest <- read.table("novel_II_ex.bed", header=F, stringsAsFactors=F)
novel_III_5 <- read.table("novel_III_5.bed", header=F, stringsAsFactors=F)
novel_III_3 <- read.table("novel_III_3.bed", header=F, stringsAsFactors=F)
novel_III_rest <- read.table("novel_III_ex.bed", header=F, stringsAsFactors=F)
intergenic_5 <- read.table("intergenic_5.bed", header=F, stringsAsFactors=F)
intergenic_3 <- read.table("intergenic_3.bed", header=F, stringsAsFactors=F)
intergenic_rest <- read.table("intergenic_ex.bed", header=F, stringsAsFactors=F)
#This gives you number of occurances when 
novel_I_5_sum <- sum(novel_I_5[["V13"]]) #13969
novel_I_3_sum <- sum(novel_I_3[["V13"]]) #13969 
novel_II_5_sum <- sum(novel_II_5[["V13"]]) #3895 
novel_II_3_sum <- sum(novel_II_3[["V13"]]) #3548
novel_III_5_sum <- sum(novel_III_5[["V13"]]) #1345
novel_III_3_sum <- sum(novel_III_3[["V13"]]) #1433
intergenic_5_sum <- sum(intergenic_5[["V13"]]) #3653
intergenic_3_sum <- sum(intergenic_3[["V13"]]) #3411
#Now to look at number of lncRNA falling in this area
novel_I_5_in <- subset(novel_I_5, (V13 > 0))#4792
novel_I_3_in <- subset(novel_I_3, (V13 > 0))#136
novel_I_upAnddown <- rbind(novel_I_5_in,novel_I_3_in)
novel_I_upAnddown <- novel_I_upAnddown[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")] 
rownames(novel_I_upAnddown) <- c()
novel_II_5_in <- subset(novel_II_5, (V13 > 0))#2056
novel_II_3_in <- subset(novel_II_3, (V13 > 0))#53
novel_II_upAnddown <- rbind(novel_II_5_in,novel_II_3_in)
novel_II_upAnddown <- novel_II_upAnddown[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")] 
rownames(novel_II_upAnddown) <- c()
novel_III_5_in <- subset(novel_III_5, (V13 > 0))#866
novel_III_3_in <- subset(novel_III_3, (V13 > 0))#22
novel_III_upAnddown <- rbind(novel_III_5_in,novel_III_3_in)
novel_III_upAnddown <- novel_III_upAnddown[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")] 
rownames(novel_III_upAnddown) <- c()
intergenic_5_in <- subset(intergenic_5, (V13 > 0))#2015
intergenic_3_in <- subset(intergenic_3, (V13 > 0))#66
intergenic_upAnddown <- rbind(intergenic_5_in,intergenic_3_in)
intergenic_upAnddown <- intergenic_upAnddown[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")] 
rownames(intergenic_upAnddown) <- c()

#remove 3' upstream and 5' upstream lncRNA
require(dplyr)
novel_I_bed_f3 <- anti_join(novel_I_bed,novel_I_upAnddown, by="V4")
novel_II_bed_f3 <- anti_join(novel_II_bed,novel_II_upAnddown, by="V4")
novel_III_bed_f3 <- anti_join(novel_III_bed,novel_III_upAnddown, by="V4")
intergenic_bed_f3 <- anti_join(intergenic_bed,intergenic_upAnddown, by="V4")
#order chr properly
novel_I_bed_f3 <- novel_I_bed_f3[with(novel_I_bed_f3, order(V1, V2)), ]
novel_II_bed_f3 <- novel_II_bed_f3[with(novel_II_bed_f3, order(V1, V2)), ]
novel_III_bed_f3 <- novel_III_bed_f3[with(novel_III_bed_f3, order(V1, V2)), ]
intergenic_bed_f3 <- intergenic_bed_f3[with(intergenic_bed_f3, order(V1, V2)), ]
#write ne bed file
setwd("~/Desktop/lncRNA")
write.table(novel_I_bed_f3, "novel_I_f3.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_bed_f3, "novel_II_f3.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_bed_f3, "novel_III_f3.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_bed_f3, "intergenic_f3.bed", row.names=F, col.names=F, quote=F, sep = "\t")

###Looking at BLASTp results
##looking at blastp.outfmt6 tables
novel_I_blastp <- read.table("novel_I_blastp.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
novel_II_blastp <- read.table("novel_II_blastp.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
novel_III_blastp <- read.table("novel_III_blastp.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
intergenic_blastp <- read.table("intergenic_blastp.outfmt6", header=F, stringsAsFactors=F,sep = "\t")
#parse out column V1 to get gene IDs
require(tidyr)
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

#write blastp tables into comparable format with hmmersearch results
setwd("~/Desktop/lncRNA")
write.table(novel_I_blastp_bed, "novel_I_blastp.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_blastp_bed, "novel_II_blastp.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_blastp_bed, "novel_III_blastp.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_blastp_bed, "intergenic_blastp.bed", row.names=F, col.names=F, quote=F, sep = "\t")

##look at HMMER PFAM results
setwd("~/Desktop/lncRNA")
novel_I_pfam <- read.table("novel_I_pfam_new.tblout", header=F, stringsAsFactors=F)
novel_II_pfam <- read.table("novel_II_pfam_new.tblout", header=F, stringsAsFactors=F)
novel_III_pfam <- read.table("novel_III_pfam_new.tblout", header=F, stringsAsFactors=F)
intergenic_pfam <- read.table("intergenic_pfam_new.tblout", header=F, stringsAsFactors=F)
###processing hmmsearch .tblout output to be compatible with blastp
require(tidyr)
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

setwd("~/Desktop/lncRNA")
write.table(novel_I_pfam_sub, "novel_I_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_pfam_sub, "novel_II_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_pfam_sub, "novel_III_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_pfam_sub, "intergenic_pfam_sub.bed", row.names=F, col.names=F, quote=F, sep = "\t")
# order of table: chr,starts,stop,V3-Name,length

#format tables for bedtools intersect
names(novel_I_blastp_bed)<-names(novel_I_pfam_sub)
names(novel_II_blastp_bed)<-names(novel_II_pfam_sub)
names(novel_III_blastp_bed)<-names(novel_III_pfam_sub)
names(intergenic_blastp_bed)<-names(intergenic_pfam_sub)

novel_I_bed_trunc <-novel_I_bed_f3[ ,c("V1","V2","V3","V4")]
novel_II_bed_trunc <-novel_II_bed_f3[ ,c("V1","V2","V3","V4")]
novel_III_bed_trunc <-novel_III_bed_f3[ ,c("V1","V2","V3","V4")]
intergenic_bed_trunc <- intergenic_bed_f3[ ,c("V1","V2","V3","V4")]

novel_I_pfam_sub_trunc <-novel_I_pfam_sub[ ,c("chr","start","stop","V3")]
novel_II_pfam_sub_trunc <-novel_II_pfam_sub[ ,c("chr","start","stop","V3")]
novel_III_pfam_sub_trunc <-novel_III_pfam_sub[ ,c("chr","start","stop","V3")]
intergenic_pfam_sub_trunc <-intergenic_pfam_sub[ ,c("chr","start","stop","V3")]

novel_I_blastp_bed_trunc <-novel_I_blastp_bed[ ,c("chr","start","stop","V3")]
novel_II_blastp_bed_trunc <-novel_II_blastp_bed[ ,c("chr","start","stop","V3")]
novel_III_blastp_bed_trunc <-novel_III_blastp_bed[ ,c("chr","start","stop","V3")]
intergenic_blastp_bed_trunc <-intergenic_blastp_bed[ ,c("chr","start","stop","V3")]

#order chr properly
novel_I_bed_trunc <- novel_I_bed_trunc[with(novel_I_bed_trunc, order(V1, V2)), ]
novel_II_bed_trunc <- novel_II_bed_trunc[with(novel_II_bed_trunc, order(V1, V2)), ]
novel_III_bed_trunc <- novel_III_bed_trunc[with(novel_III_bed_trunc, order(V1, V2)), ]
intergenic_bed_trunc <- intergenic_bed_trunc[with(intergenic_bed_trunc, order(V1, V2)), ]

novel_I_pfam_sub_trunc <- novel_I_pfam_sub_trunc[with(novel_I_pfam_sub_trunc, order(chr, start)), ]
novel_II_pfam_sub_trunc <- novel_II_pfam_sub_trunc[with(novel_II_pfam_sub_trunc, order(chr, start)), ]
novel_III_pfam_sub_trunc <- novel_III_pfam_sub_trunc[with(novel_III_pfam_sub_trunc, order(chr, start)), ]
intergenic_pfam_sub_trunc <- intergenic_pfam_sub_trunc[with(intergenic_pfam_sub_trunc, order(chr, start)), ]

novel_I_blastp_bed_trunc <- novel_I_blastp_bed_trunc[with(novel_I_blastp_bed_trunc, order(chr, start)), ]
novel_II_blastp_bed_trunc <- novel_II_blastp_bed_trunc[with(novel_II_blastp_bed_trunc, order(chr, start)), ]
novel_III_blastp_bed_trunc <- novel_III_blastp_bed_trunc[with(novel_III_blastp_bed_trunc, order(chr, start)), ]
intergenic_blastp_bed_trunc <- intergenic_blastp_bed_trunc[with(intergenic_blastp_bed_trunc, order(chr, start)), ]


write.table(novel_I_bed_trunc, "novel_I_trunc.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_bed_trunc, "novel_II_trunc.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_bed_trunc, "novel_III_trunc.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_bed_trunc, "intergenic_trunc.bed", row.names=F, col.names=F, quote=F, sep = "\t")
#Hmmer tables
write.table(novel_I_pfam_sub_trunc, "novel_I_H_trunc.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_pfam_sub_trunc, "novel_II_H_trunc.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_pfam_sub_trunc, "novel_III_H_trunc.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_pfam_sub_trunc, "intergenic_H_trunc.bed", row.names=F, col.names=F, quote=F, sep = "\t")
#Blastp tables
write.table(novel_I_blastp_bed_trunc, "novel_I_B_trunc.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_blastp_bed_trunc, "novel_II_B_trunc.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_blastp_bed_trunc, "novel_III_B_trunc.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_blastp_bed_trunc, "intergenic_B_trunc.bed", row.names=F, col.names=F, quote=F, sep = "\t")

#format tables for bedtools intersect
novel_I_bed_trunc <-novel_I_bed_f3[ ,c("V1","V2","V3","V4")]
novel_II_bed_trunc <-novel_II_bed_f3[ ,c("V1","V2","V3","V4")]
novel_III_bed_trunc <-novel_III_bed_f3[ ,c("V1","V2","V3","V4")]
intergenic_bed_trunc <- intergenic_bed_f3[ ,c("V1","V2","V3","V4")]
##anti-join to get non PCG
names(novel_I_blastp_bed)<-names(novel_I_bed_trunc)
names(novel_II_blastp_bed)<-names(novel_II_bed_trunc)
names(novel_III_blastp_bed)<-names(novel_III_bed_trunc)
names(intergenic_blastp_bed)<-names(intergenic_bed_trunc)

names(novel_I_pfam_sub_trunc)<-names(novel_I_bed_trunc)
names(novel_II_pfam_sub_trunc)<-names(novel_II_bed_trunc)
names(novel_III_pfam_sub_trunc)<-names(novel_III_bed_trunc)
names(intergenic_pfam_sub_trunc)<-names(intergenic_bed_trunc)
##filtering out PCGs 
novel_I_P_bed <- merge(novel_I_bed_trunc,novel_I_pfam_sub_trunc, by.x=c("V1","V2","V3"),by.y=c("V1","V2","V3"))
novel_I_noP_bed <- anti_join(novel_I_bed_trunc,novel_I_P_bed, by.x="V4",by.y="V4.x")

novel_II_P_bed <- merge(novel_II_bed_trunc,novel_II_pfam_sub_trunc, by.x=c("V1","V2","V3"),by.y=c("V1","V2","V3"))
novel_II_noP_bed <- anti_join(novel_II_bed_trunc,novel_II_P_bed, by.x="V4",by.y="V4.x")

novel_III_P_bed <- merge(novel_III_bed_trunc,novel_III_pfam_sub_trunc, by.x=c("V1","V2","V3"),by.y=c("V1","V2","V3"))
novel_III_noP_bed <- anti_join(novel_III_bed_trunc,novel_III_P_bed, by.x="V4",by.y="V4.x")

intergenic_P_bed <- merge(intergenic_bed_trunc,intergenic_pfam_sub_trunc, by.x=c("V1","V2","V3"),by.y=c("V1","V2","V3"))
intergenic_noP_bed <- anti_join(intergenic_bed_trunc,intergenic_P_bed, by.x="V4",by.y="V4.x")


write.table(novel_I_P_bed, "novel_I_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_I_noP_bed, "novel_I_noP.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_P_bed, "novel_II_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_noP_bed, "novel_II_noP.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_P_bed, "novel_III_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_noP_bed, "novel_III_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_P_bed, "intergenic_P.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_noP_bed, "intergenic_noP.bed", row.names=F, col.names=F, quote=F, sep = "\t")

all_lncRNA <- rbind(novel_I_noP_bed,novel_II_noP_bed,novel_III_noP_bed,intergenic_noP_bed)
all_lncRNA <- all_lncRNA[with(all_lncRNA, order(V1, V2)), ]
lncRNA_all_Cat <-rbind(data.frame(id="novel_I",novel_I_noP_bed),
                   data.frame(id="novel_II",novel_II_noP_bed),
                   data.frame(id="novel_III",novel_III_noP_bed),
                   data.frame(id="intergenic",intergenic_noP_bed))
lncRNA_all_Cat <- lncRNA_all_Cat[with(lncRNA_all_Cat, order(V1, V2)), ]
write.table(all_lncRNA, "lncRNA_final8247", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(lncRNA_all_Cat, "lncRNA_final8247_IDs", row.names=F, col.names=F, quote=F, sep = "\t")

novel_I_P_bed <- novel_I_P_bed[ ,c("V1","V2","V3","V4.x")]
novel_I_P_bed <- novel_I_P_bed[!duplicated(novel_I_P_bed),]
names(novel_I_P_bed)[4]<-paste("V4")

novel_II_P_bed <- novel_II_P_bed[ ,c("V1","V2","V3","V4.x")]
novel_II_P_bed <- novel_II_P_bed[!duplicated(novel_II_P_bed),]
names(novel_II_P_bed)[4]<-paste("V4")

novel_III_P_bed <- novel_III_P_bed[ ,c("V1","V2","V3","V4.x")]
novel_III_P_bed <- novel_III_P_bed[!duplicated(novel_III_P_bed),]
names(novel_III_P_bed)[4]<-paste("V4")

intergenic_P_bed <- intergenic_P_bed[ ,c("V1","V2","V3","V4.x")]
intergenic_P_bed <- intergenic_P_bed[!duplicated(intergenic_P_bed),]
names(intergenic_P_bed)[4]<-paste("V4")


all_ID <-rbind(data.frame(id="novel_I_lncRNA",novel_I_noP_bed),
               data.frame(id="novel_I_genes",novel_I_P_bed),
               data.frame(id="novel_II_lncRNA",novel_II_noP_bed),
               data.frame(id="novel_II_genes",novel_II_P_bed),
               data.frame(id="novel_III_lncRNA",novel_III_noP_bed),
               data.frame(id="novel_III_genes",novel_III_P_bed),
               data.frame(id="intergenic_lncRNA",intergenic_noP_bed),
               data.frame(id="intergenic_genes",intergenic_P_bed))

write.table(all_ID, "all_cats_PandnoP", row.names=F, col.names=F, quote=F, sep = "\t")


