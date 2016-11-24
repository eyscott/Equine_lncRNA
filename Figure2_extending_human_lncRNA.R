##extending boundaries on human lncRNA bed file
human_lncRNA_bed <- read.table("hg19.lncRNAs.bed", header=F, stringsAsFactors=F)
human_lncRNA_bed["up1"] <- (human_lncRNA_bed[ ,2]-1000)
human_lncRNA_bed["down1"]<- (human_lncRNA_bed[ ,3] + 1000)
human_lncRNA_bed_1 <- human_lncRNA_bed[ ,c("V1","up1","down1","V4","V5","V6","V7","V8","V9","V10","V11","V12")]

human_lncRNA_bed["up2"] <- (human_lncRNA_bed[ ,2]-2000)
human_lncRNA_bed["down2"]<- (human_lncRNA_bed[ ,3] + 2000)
human_lncRNA_bed_2 <- human_lncRNA_bed[ ,c("V1","up2","down2","V4","V5","V6","V7","V8","V9","V10","V11","V12")]

human_lncRNA_bed["up3"] <- (human_lncRNA_bed[ ,2]-3000)
human_lncRNA_bed["down3"]<- (human_lncRNA_bed[ ,3] + 3000)
human_lncRNA_bed_3 <- human_lncRNA_bed[ ,c("V1","up3","down3","V4","V5","V6","V7","V8","V9","V10","V11","V12")]

human_lncRNA_bed["up4"] <- (human_lncRNA_bed[ ,2]-4000)
human_lncRNA_bed["down4"]<- (human_lncRNA_bed[ ,3] + 4000)
human_lncRNA_bed_4 <- human_lncRNA_bed[ ,c("V1","up4","down4","V4","V5","V6","V7","V8","V9","V10","V11","V12")]

human_lncRNA_bed["up5"] <- (human_lncRNA_bed[ ,2]-5000)
human_lncRNA_bed["down5"]<- (human_lncRNA_bed[ ,3] + 5000)
human_lncRNA_bed_5 <- human_lncRNA_bed[ ,c("V1","up5","down5","V4","V5","V6","V7","V8","V9","V10","V11","V12")]

human_lncRNA_bed["up6"] <- (human_lncRNA_bed[ ,2]-6000)
human_lncRNA_bed["down6"]<- (human_lncRNA_bed[ ,3] + 6000)
human_lncRNA_bed_6 <- human_lncRNA_bed[ ,c("V1","up6","down6","V4","V5","V6","V7","V8","V9","V10","V11","V12")]

human_lncRNA_bed["up7"] <- (human_lncRNA_bed[ ,2]-7000)
human_lncRNA_bed["down7"]<- (human_lncRNA_bed[ ,3] + 7000)
human_lncRNA_bed_7 <- human_lncRNA_bed[ ,c("V1","up7","down7","V4","V5","V6","V7","V8","V9","V10","V11","V12")]

human_lncRNA_bed["up8"] <- (human_lncRNA_bed[ ,2]-8000)
human_lncRNA_bed["down8"]<- (human_lncRNA_bed[ ,3] + 8000)
human_lncRNA_bed_8 <- human_lncRNA_bed[ ,c("V1","up8","down8","V4","V5","V6","V7","V8","V9","V10","V11","V12")]

human_lncRNA_bed["up9"] <- (human_lncRNA_bed[ ,2]-9000)
human_lncRNA_bed["down9"]<- (human_lncRNA_bed[ ,3] + 9000)
human_lncRNA_bed_9 <- human_lncRNA_bed[ ,c("V1","up9","down9","V4","V5","V6","V7","V8","V9","V10","V11","V12")]

human_lncRNA_bed["up10"] <- (human_lncRNA_bed[ ,2]-10000)
human_lncRNA_bed["down10"]<- (human_lncRNA_bed[ ,3] + 10000)
human_lncRNA_bed_10 <- human_lncRNA_bed[ ,c("V1","up10","down10","V4","V5","V6","V7","V8","V9","V10","V11","V12")]


#order chr properly
human_lncRNA_bed <- human_lncRNA_bed[ ,c(1:12)]
human_lncRNA_bed <- human_lncRNA_bed[with(human_lncRNA_bed, order(V1, up1)), ]
human_lncRNA_bed_1 <- human_lncRNA_bed_1[with(human_lncRNA_bed_1, order(V1, up1)), ]
human_lncRNA_bed_2 <- human_lncRNA_bed_2[with(human_lncRNA_bed_2, order(V1, up2)), ]
human_lncRNA_bed_3 <- human_lncRNA_bed_3[with(human_lncRNA_bed_3, order(V1, up3)), ]
human_lncRNA_bed_4 <- human_lncRNA_bed_4[with(human_lncRNA_bed_4, order(V1, up4)), ]
human_lncRNA_bed_5 <- human_lncRNA_bed_5[with(human_lncRNA_bed_5, order(V1, up5)), ]
human_lncRNA_bed_6 <- human_lncRNA_bed_6[with(human_lncRNA_bed_6, order(V1, up6)), ]
human_lncRNA_bed_7 <- human_lncRNA_bed_7[with(human_lncRNA_bed_7, order(V1, up7)), ]
human_lncRNA_bed_8 <- human_lncRNA_bed_8[with(human_lncRNA_bed_8, order(V1, up8)), ]
human_lncRNA_bed_9 <- human_lncRNA_bed_9[with(human_lncRNA_bed_9, order(V1, up9)), ]
human_lncRNA_bed_10 <- human_lncRNA_bed_10[with(human_lncRNA_bed_10, order(V1, up10)), ]


write.table(human_lncRNA_bed, "hg19.lncRNAs.bed", row.names=F, col.names=F, quote=F,, sep = "\t")
write.table(human_lncRNA_bed_1, "hg19.lncRNAs_1.bed", row.names=F, col.names=F, quote=F,, sep = "\t")
write.table(human_lncRNA_bed_2, "hg19.lncRNAs_2.bed", row.names=F, col.names=F, quote=F,, sep = "\t")
write.table(human_lncRNA_bed_3, "hg19.lncRNAs_3.bed", row.names=F, col.names=F, quote=F,, sep = "\t")
write.table(human_lncRNA_bed_4, "hg19.lncRNAs_4.bed", row.names=F, col.names=F, quote=F,, sep = "\t")
write.table(human_lncRNA_bed_5, "hg19.lncRNAs_5.bed", row.names=F, col.names=F, quote=F,, sep = "\t")
write.table(human_lncRNA_bed_6, "hg19.lncRNAs_6.bed", row.names=F, col.names=F, quote=F,, sep = "\t")
write.table(human_lncRNA_bed_7, "hg19.lncRNAs_7.bed", row.names=F, col.names=F, quote=F,, sep = "\t")
write.table(human_lncRNA_bed_8, "hg19.lncRNAs_8.bed", row.names=F, col.names=F, quote=F,, sep = "\t")
write.table(human_lncRNA_bed_9, "hg19.lncRNAs_9.bed", row.names=F, col.names=F, quote=F,, sep = "\t")
write.table(human_lncRNA_bed_10, "hg19.lncRNAs_10.bed", row.names=F, col.names=F, quote=F,, sep = "\t")


#import and order chr properly for mapped horse lncRNA (to human)
horse_to_human_bed <- read.table("all_lncRNA_final_hg19_mapped.bed", header=F, stringsAsFactors=F) 
horse_to_human_bed <- horse_to_human_bed[with(horse_to_human_bed, order(V1, V2)), ]
write.table(horse_to_human_bed, "mapped_lncRNA_horseThuman.bed", row.names=F, col.names=F, quote=F,, sep = "\t")

#After running bedtools -intersect with al extended versions of human lncRNA
##read in bedtools -intersect files
setwd("~/Desktop/lncRNA")
human_bed <- read.table("lncRNA_human.bed", header=F, stringsAsFactors=F)
human_bed_1 <- read.table("lncRNA_human_1.bed", header=F, stringsAsFactors=F)
human_bed_2 <- read.table("lncRNA_human_2.bed", header=F, stringsAsFactors=F)
human_bed_3 <- read.table("lncRNA_human_3.bed", header=F, stringsAsFactors=F)
human_bed_4 <- read.table("lncRNA_human_4.bed", header=F, stringsAsFactors=F)
human_bed_5 <- read.table("lncRNA_human_5.bed", header=F, stringsAsFactors=F)
human_bed_6 <- read.table("lncRNA_human_6.bed", header=F, stringsAsFactors=F)
human_bed_7 <- read.table("lncRNA_human_7.bed", header=F, stringsAsFactors=F)
human_bed_8 <- read.table("lncRNA_human_8.bed", header=F, stringsAsFactors=F)
human_bed_9 <- read.table("lncRNA_human_9.bed", header=F, stringsAsFactors=F)
human_bed_10 <- read.table("lncRNA_human_10.bed", header=F, stringsAsFactors=F)

human_bed_in <- subset(human_bed, (V13 > 0))#2448
human_bed_1_in <- subset(human_bed_1, (V13 > 0))#2623
human_bed_2_in <- subset(human_bed_2, (V13 > 0))#2752
human_bed_3_in <- subset(human_bed_3, (V13 > 0))#2865
human_bed_4_in <- subset(human_bed_4, (V13 > 0))#2969
human_bed_5_in <- subset(human_bed_5, (V13 > 0))#3078
human_bed_6_in <- subset(human_bed_6, (V13 > 0))#3168
human_bed_7_in <- subset(human_bed_7, (V13 > 0))#3254
human_bed_8_in <- subset(human_bed_8, (V13 > 0))#3346
human_bed_9_in <- subset(human_bed_9, (V13 > 0))#3425
human_bed_10_in <- subset(human_bed_10, (V13 > 0))#3506


d <- c(10,9,8,7,6,5,4,3,2,1,0,1,2,3,4,5,6,7,8,9,10)
N_lncRNA <- c(3506,3425,3346,3254,3168,3078,2969,2865,2752,2623,2448,2623,2752,2865,2969,3078,3168,3254,3346,3425,3506)
dist_ncRNA <- data.frame(d,N_lncRNA)

#plotting the distances from 9402 lncRNA mapped to human 
library(ggplot2)
ggplot(dist_ncRNA,aes(d,N_lncRNA)) +
  geom_line() + ylab("Number of candidate lncRNA") +
  xlab("distance (kb) to an annotated human lncRNA")

