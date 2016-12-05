##transcripts to 5' and 3'...reading the bedtools intersect files
setwd("~/Desktop/lncRNA")
novel_I_5 <- read.table("novel_I_5.bed", header=F, stringsAsFactors=F)
novel_I_3 <- read.table("novel_I_3.bed", header=F, stringsAsFactors=F)
novel_II_5 <- read.table("novel_II_5.bed", header=F, stringsAsFactors=F)
novel_II_3 <- read.table("novel_II_3.bed", header=F, stringsAsFactors=F)
novel_III_5 <- read.table("novel_III_5.bed", header=F, stringsAsFactors=F)
novel_III_3 <- read.table("novel_III_3.bed", header=F, stringsAsFactors=F)
intergenic_5 <- read.table("intergenic_5.bed", header=F, stringsAsFactors=F)
intergenic_3 <- read.table("intergenic_3.bed", header=F, stringsAsFactors=F)
known_5 <- read.table("known_5.bed", header=F, stringsAsFactors=F)
known_3 <- read.table("known_3.bed", header=F, stringsAsFactors=F)

#look at number of lncRNA falling in this area
novel_I_5_in <- subset(novel_I_5, (V13 > 0))#4610
novel_I_3_in <- subset(novel_I_3, (V13 > 0))#4662
novel_I_upAnddown <- rbind(novel_I_5_in,novel_I_3_in)
novel_I_upAnddown <- novel_I_upAnddown[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")] 
rownames(novel_I_upAnddown) <- c()
#novel_I_upAnddown are the F3 rejects for novel I
write.table(novel_I_upAnddown, "F4_novel_I", row.names=F, col.names=F, quote=F, sep = "\t")
novel_II_5_in <- subset(novel_II_5, (V13 > 0))#1487
novel_II_3_in <- subset(novel_II_3, (V13 > 0))#1418
novel_II_upAnddown <- rbind(novel_II_5_in,novel_II_3_in)
novel_II_upAnddown <- novel_II_upAnddown[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")] 
rownames(novel_II_upAnddown) <- c()
#novel_II_upAnddown are the F3 rejects for novel II
write.table(novel_II_upAnddown, "F4_novel_II", row.names=F, col.names=F, quote=F, sep = "\t")
novel_III_5_in <- subset(novel_III_5, (V13 > 0))#268
novel_III_3_in <- subset(novel_III_3, (V13 > 0))#276
novel_III_upAnddown <- rbind(novel_III_5_in,novel_III_3_in)
novel_III_upAnddown <- novel_III_upAnddown[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")] 
rownames(novel_III_upAnddown) <- c()
#novel_III_upAnddown are the F3 rejects for novel III
write.table(novel_III_upAnddown, "F4_novel_III", row.names=F, col.names=F, quote=F, sep = "\t")
intergenic_5_in <- subset(intergenic_5, (V13 > 0))#2392
intergenic_3_in <- subset(intergenic_3, (V13 > 0))#2341
intergenic_upAnddown <- rbind(intergenic_5_in,intergenic_3_in)
intergenic_upAnddown <- intergenic_upAnddown[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")] 
rownames(intergenic_upAnddown) <- c()
d <- intergenic_upAnddown[!duplicated(intergenic_upAnddown),]#4205, therefore 528 in both 5' and 3' up- and downstream sequences
#intergenic_upAnddown are the F3 rejects for intergenic
write.table(intergenic_upAnddown, "F4_intergenic", row.names=F, col.names=F, quote=F, sep = "\t")
known_5_in <- subset(known_5, (V13 > 0))#2392
known_3_in <- subset(known_3, (V13 > 0))#2341
known_upAnddown <- rbind(known_5_in,known_3_in)
known_upAnddown <- known_upAnddown[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")] 
rownames(known_upAnddown) <- c()
#known_upAnddown are the F3 rejects for intergenic
write.table(known_upAnddown, "F4_known", row.names=F, col.names=F, quote=F, sep = "\t")



#remove 3' upstream and 5' upstream lncRNA
setwd("~/Desktop/lncRNA")
novel_I_bed <-read.table("novel_I_f3.bed",header=F)
novel_II_bed <-read.table("novel_II_f3.bed",header=F)
novel_III_bed <-read.table("novel_III_f3.bed",header=F)
intergenic_bed <-read.table("intergenic_f3.bed",header=F)
known_bed <-read.table("known_lncRNA_f3.bed",header=F)

require(dplyr)
novel_I_bed_f4 <- anti_join(novel_I_bed,novel_I_upAnddown, by="V4")
novel_II_bed_f4 <- anti_join(novel_II_bed,novel_II_upAnddown, by="V4")
novel_III_bed_f4 <- anti_join(novel_III_bed,novel_III_upAnddown, by="V4")
intergenic_bed_f4 <- anti_join(intergenic_bed,intergenic_upAnddown, by="V4")
known_bed_f4 <- anti_join(known_bed,known_upAnddown, by="V4")

#order chr properly
novel_I_bed_f4 <- novel_I_bed_f4[with(novel_I_bed_f4, order(V1, V2)), ]
novel_II_bed_f4 <- novel_II_bed_f4[with(novel_II_bed_f4, order(V1, V2)), ]
novel_III_bed_f4 <- novel_III_bed_f4[with(novel_III_bed_f4, order(V1, V2)), ]
intergenic_bed_f4 <- intergenic_bed_f4[with(intergenic_bed_f4, order(V1, V2)), ]
known_bed_f4 <- known_bed_f4[with(known_bed_f4, order(V1, V2)), ]

#write the bed file
setwd("~/Desktop/lncRNA")
write.table(novel_I_bed_f4, "novel_I_f4.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_bed_f4, "novel_II_f4.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_bed_f4, "novel_III_f4.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_bed_f4, "intergenic_f4.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(kown_bed_f4, "intergenic_f4.bed", row.names=F, col.names=F, quote=F, sep = "\t")
