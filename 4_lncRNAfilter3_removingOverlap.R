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
novel_I_5_in <- subset(novel_I_5, (V13 > 0))#4610
novel_I_3_in <- subset(novel_I_3, (V13 > 0))#4662
novel_I_upAnddown <- rbind(novel_I_5_in,novel_I_3_in)
novel_I_upAnddown <- novel_I_upAnddown[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")] 
rownames(novel_I_upAnddown) <- c()
#novel_I_upAnddown are the F3 rejects for novel I
write.table(novel_I_upAnddown, "F3_novel_I", row.names=F, col.names=F, quote=F, sep = "\t")
novel_II_5_in <- subset(novel_II_5, (V13 > 0))#1487
novel_II_3_in <- subset(novel_II_3, (V13 > 0))#1418
novel_II_upAnddown <- rbind(novel_II_5_in,novel_II_3_in)
novel_II_upAnddown <- novel_II_upAnddown[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")] 
rownames(novel_II_upAnddown) <- c()
#novel_II_upAnddown are the F3 rejects for novel II
write.table(novel_II_upAnddown, "F3_novel_II", row.names=F, col.names=F, quote=F, sep = "\t")
novel_III_5_in <- subset(novel_III_5, (V13 > 0))#268
novel_III_3_in <- subset(novel_III_3, (V13 > 0))#276
novel_III_upAnddown <- rbind(novel_III_5_in,novel_III_3_in)
novel_III_upAnddown <- novel_III_upAnddown[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")] 
rownames(novel_III_upAnddown) <- c()
#novel_III_upAnddown are the F3 rejects for novel III
write.table(novel_III_upAnddown, "F3_novel_III", row.names=F, col.names=F, quote=F, sep = "\t")
intergenic_5_in <- subset(intergenic_5, (V13 > 0))#2392
intergenic_3_in <- subset(intergenic_3, (V13 > 0))#2341
intergenic_upAnddown <- rbind(intergenic_5_in,intergenic_3_in)
intergenic_upAnddown <- intergenic_upAnddown[ ,c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")] 
rownames(intergenic_upAnddown) <- c()
d <- intergenic_upAnddown[!duplicated(intergenic_upAnddown),]#2824, therefore 19 in both 5' and 3' U
#intergenic_upAnddown are the F3 rejects for intergenic
write.table(intergenic_upAnddown, "F3_intergenic", row.names=F, col.names=F, quote=F, sep = "\t")
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
#write the bed file
setwd("~/Desktop/lncRNA")
write.table(novel_I_bed_f3, "novel_I_f3_new.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_bed_f3, "novel_II_f3_new.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_bed_f3, "novel_III_f3_new.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_bed_f3, "intergenic_f3_new.bed", row.names=F, col.names=F, quote=F, sep = "\t")

