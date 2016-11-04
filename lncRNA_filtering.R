##lncRNA filtering
##load in initial backmapping TPM data, added transcriptName label to TCONS IDs
setwd("~/Desktop/lncRNA")
initial_data <- read.table("dataSummary", header=T, stringsAsFactors=F)

#filter out transcripts <200 nt and with TPM < 0.1
f1 <- subset(initial_data, c(calcTPM > 0.1 & length > 200))

##remove any transcripts already listed in our reference transcriptome
#use the unfiltered bed to get genomic coordinates
unfiltered_bed <- read.table("unfiltered_Alltissues_Assembly.bed", header=F, stringsAsFactors=F)
m <- merge(f1, unfiltered_bed, by.x="transcriptName",by.y="V4" )
lncRNA_bed <- m[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(lncRNA_bed)[4]<-paste("V4")
#Use the refined bed for removal of genes overlapping genes
refined_bed <- read.table("refined.bed", header=F,stringsAsFactors=F )
#Use the merged bed for removal of overlapping genes
merged_bed <- read.table("mergedTrans.bed", header=F,stringsAsFactors=F )
#find the unique transcript ID's not present in the refined.bed, aka non-annotated transcripts
require(dplyr) 
f2 <- anti_join(lncRNA_bed,refined_bed, by="V4")
f2_merged <- anti_join(lncRNA_bed,merged_bed, by="V4")
#brings our db down to 64668 trancripts, with merged=67561

#minus 1 kb and add 1 kb to end of each gene in refined.bed
refined_bed["minus"] <- (refined_bed[ ,2]-1000)
refined_bed["add"]<- (refined_bed[ ,3] + 1000)
refined_bed_extended <- refined_bed[ ,c("V1","minus","add","V4","V5","V6","V7","V8","V9","V10","V11","V12")]
names(refined_bed_extended)[2]<-paste("V2")
names(refined_bed_extended)[3]<-paste("V3")
#now with merged bed
merged_bed["minus"] <- (merged_bed[ ,2]-1000)
merged_bed["add"]<- (merged_bed[ ,3] + 1000)
merged_bed_extended <- merged_bed[ ,c("V1","minus","add","V4","V5","V6","V7","V8","V9","V10","V11","V12")]
names(merged_bed_extended)[2]<-paste("V2")
names(merged_bed_extended)[3]<-paste("V3")
#this is to allow for removal of any transcripts within 1 kb of an annotated gene

#convert all negative numbers into 0
refined_bed_extended$V2[refined_bed_extended$V2 < 0] <- 0
merged_bed_extended$V2[merged_bed_extended$V2 < 0] <- 0

#sort bed files by chr then by genomic coordinates
lncRNA_bed_f2 <- f2[with(f2, order(V1, V2)), ]
lncRNA_bed_f2_merged <- f2_merged[with(f2_merged, order(V1, V2)), ]
refined_bed_extended <- refined_bed_extended[with(refined_bed_extended, order(V1, V2)), ]
merged_bed_extended <- merged_bed_extended[with(merged_bed_extended, order(V1, V2)), ]

#write tables making sure that column and row names are not printed as well...
write.table(refined_bed_extended, "refined_extended.bed", row.names=F, col.names=F)
write.table(lncRNA_bed_f2, "lncRNA.bed", row.names=F, col.names=F)
write.table(lncRNA_bed_f2_merged, "lncRNA_merged.bed", row.names=F, col.names=F)
write.table(merged_bed_extended, "merged_extended.bed", row.names=F, col.names=F)

##use bedtool intersect function done on other server
tr '\r' '\n' < lncRNA.bed > lncRNA_unix.bed 
tr '\r' '\n' < merged_extended.bed > merged_extended_unix.bed
tr '\r' '\n' < lncRNA_merged.bed > lncRNA_merged_unix.bed 
tr '\r' '\n' < refined_extended.bed > refined_extended_unix.bed
tr '\r' '\n' < mergedTrans.bed > mergedTrans_unix.bed

bedtools intersect  -c -v -sorted -a lncRNA_unix.bed -b refined_extended_unix.bed > cv.bed
bedtools intersect  -v -sorted -a lncRNA_unix.bed -b refined_extended_unix.bed > v.bed
bedtools intersect  -c -v -sorted -a lncRNA_merged_unix.bed -b merged_extended.bed > merged_cv.bed
bedtools intersect  -c -sorted -a lncRNA_unix.bed -b refined_extended_unix.bed > c.bed

##result:
cv_bed <- read.table("cv.bed", header=F, stringsAsFactors=F)
##getting expression values for lncRNA transcripts
lncRNA_bed<- read.table("lncRNA.bed", header=F, stringsAsFactors=F)
lncRNA_exp<- merge(lncRNA_bed, f1, by.x="V4", by.y="transcriptName")
mean_TPM<-mean(lncRNA_exp[["calcTPM"]])

##after bedtools....
setwd("~/Desktop/lncRNA/NEW")
v_bed <- read.table("v.bed", header=F, stringsAsFactors=F)
c_bed <- read.table("c.bed", header=F, stringsAsFactors=F)
mergedv_bed <- read.table("merged_v.bed", header=F, stringsAsFactors=F)

library(ggplot2)
keeps <- c("V1", "V2", "V3", "V4")
v_bed_keeps <- v_bed[keeps]
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chr31", "chrUn", "chrX")
chrs_N <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "Un", "X")    
chromSizes<-read.table("equCab2.chrom.sizes.txt", header=T, col.names=c("chr","size"),stringsAsFactors=FALSE) 

m <- ggplot(data=v_bed_keeps, aes(V1)) + geom_bar(stat = "count") + scale_x_discrete(limits = chrs, labels = chrs_N) + ylab("transcript count") + xlab("Chr") +
  geom_line(data=subset(chromSizes,chr %in% chrs), aes(x=chr, y= size / 50000, group=1), colour="blue")

m

mv_bed_keeps <- mergedv_bed[keeps]
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chr24", "chr25", "chr26", "chr27", "chr28", "chr29", "chr30", "chr31", "chrUn", "chrX")
chrs_N <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "Un", "X")    
chromSizes<-read.table("equCab2.chrom.sizes.txt", header=T, col.names=c("chr","size"),stringsAsFactors=FALSE) 

n <- ggplot(data=mv_bed_keeps, aes(V1)) + geom_bar(stat = "count") + scale_x_discrete(limits = chrs, labels = chrs_N) + ylab("transcript count") + xlab("Chr") +
  geom_line(data=subset(chromSizes,chr %in% chrs), aes(x=chr, y= size / 50000, group=1), colour="blue")

n

###finding exon numbers for v_bed files
m2 <- merge(v_bed_keeps, unfiltered_bed, by.x="V4",by.y="V4" )
lncRNA_bed_bedtools <- m2[ ,c("V1.x", "V2.x",  "V3.x","V4",  "V1.y",  "V2.y",	"V3.y",	"V5",	"V6",	"V7",	"V8",	"V9",	"V10",	"V11",	"V12")] 
lncRNA_bed_bedtools2 <- lncRNA_bed_bedtools[with(lncRNA_bed_bedtools, order(V1.x, V2.x)), ]

n2 <- merge(mv_bed_keeps, unfiltered_bed, by.x="V4",by.y="V4" )
lncRNA_bed_bedtools <- n2[ ,c("V1.x", "V2.x",  "V3.x","V4",  "V1.y",	"V2.y",	"V3.y",	"V5",	"V6",	"V7",	"V8",	"V9",	"V10",	"V11",	"V12")] 
lncRNA_bed_bedtools2 <- lncRNA_bed_bedtools[with(lncRNA_bed_bedtools, order(V1.x, V2.x)), ]

write.table(lncRNA_bed_bedtools2, "lncRNA_bedtools.bed",row.names=F, col.names=F)
