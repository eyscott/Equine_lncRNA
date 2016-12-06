##lncRNA pipeline R##
setwd("~/lncRNA")
overallExpression <- read.table("inputs/dataSummary", header=T, stringsAsFactors=F)
overallExpression$transcriptName=rownames(overallExpression)

#load the different transcript categories
known_ncRNA <- read.table("RNAseqSupTrans.merge.reduced.ncRNA", header=T, stringsAsFactors=F)
novel_I <- read.table("inputs/novelAnn/sup/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.sup", header=T, stringsAsFactors=F)
novel_II <- read.table("inputs/novelAnn/unsup.cons/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.unsup.cons", header=T, stringsAsFactors=F)
novel_III<- read.table("inputs/novelAnn/unsup.uncons.ORF/RNAseqSupTrans.merge.reduced.ORF_exons.candNovel.unsup.uncons.ORF", header=T, stringsAsFactors=F)
intergenic<- read.table("inputs/intergenic/unsup.ids", header=F, stringsAsFactors=F)

#attach global gene expression values and length 
known_ncRNA <- merge(overallExpression, known_ncRNA, by.x="transcriptName",by.y="transcript.ID" )
novel_I <- merge(overallExpression, novel_I, by.x="transcriptName",by.y="transcript.ID" )
novel_II <- merge(overallExpression, novel_II, by.x="transcriptName",by.y="transcript.ID" )
novel_III <- merge(overallExpression, novel_III, by.x="transcriptName",by.y="transcript.ID" )
intergenic <- merge(overallExpression, intergenic, by.x="transcriptName",by.y="V1" )

#slim down to columns you want
keeps <- c("transcriptName","length", "calcTPM")
known_ncRNA <- unique(known_ncRNA[keeps])
novel_I <- unique(novel_I[keeps])
novel_II <- unique(novel_II[keeps])
novel_III <- unique(novel_III[keeps])
intergenic <- unique(intergenic[keeps])

#Apply initial relaxed threshold of average 0.1 TPM
require(dplyr)
#these are the subsets of transcripts passing the 0.1 TPM filter
known_ncRNA_f1 <- subset(known_ncRNA, calcTPM > 0.1)
novel_I_f1 <- subset(novel_I, calcTPM > 0.1)
novel_II_f1 <- subset(novel_II, calcTPM > 0.1)
novel_III_f1 <- subset(novel_III, calcTPM > 0.1)
intergenic_f1 <- subset(intergenic, calcTPM > 0.1)
#These are the transcripts not passing the above filter (original list- transcripts surving F1)
f1_known_rejects <- anti_join(known_ncRNA, known_ncRNA_f1, by="transcriptName")
f1_I_rejects <- anti_join(novel_I, novel_I_f1, by="transcriptName")
f1_II_rejects <- anti_join(novel_II,novel_II_f1,by="transcriptName")
f1_III_rejects <- anti_join(novel_III,novel_III_f1,by="transcriptName")
f1_intergenic_rejects <- anti_join(intergenic,intergenic_f1,by="transcriptName")

#Apply the filter 2, removing any transcripts less than 200 nt
#These are the transcripts that did pass F2
known_ncRNA_f2 <- subset(known_ncRNA_f1, length > 200)
novel_I_f2 <- subset(novel_I_f1, length > 200)
novel_II_f2 <- subset(novel_II_f1, length > 200)
novel_III_f2 <- subset(novel_III_f1, length > 200)
intergenic_f2 <- subset(intergenic_f1, length > 200)
#These are the transcripts that did not pass the F2
f2_known_rejects <- anti_join(known_ncRNA_f1, known_ncRNA_f2, by="transcriptName")
f2_I_rejects <- anti_join(novel_I_f1, novel_I_f2, by="transcriptName")
f2_II_rejects <- anti_join(novel_II_f1,novel_II_f2,by="transcriptName")
f2_III_rejects <- anti_join(novel_III_f1,novel_III_f2,by="transcriptName")
f2_intergenic_rejects <- anti_join(intergenic_f1,intergenic_f2,by="transcriptName")

##Apply more stringent criteria to single exon transcripts
#Attach exon informatin to the above inputs
#Make each into a bed file to use with bedtools
unfiltered_bed <- read.table("inputs/allTissues_BED/unfiltered_Alltissues_Assembly.bed", header=F, stringsAsFactors=F)
known_ncRNA_bed <- merge(known_ncRNA_f2, unfiltered_bed, by.x="transcriptName",by.y="V4" )
novel_I_bed <- merge(novel_I_f2, unfiltered_bed, by.x="transcriptName",by.y="V4" )
novel_II_bed <- merge(novel_II_f2, unfiltered_bed, by.x="transcriptName",by.y="V4" )
novel_III_bed <- merge(novel_III_f2, unfiltered_bed, by.x="transcriptName",by.y="V4" )
intergenic_bed <- merge(intergenic_f2, unfiltered_bed, by.x="transcriptName",by.y="V4" )

#Isolate the single exon transcripts
single_known <- subset(known_ncRNA_bed, c(V10 <2))
single_I <- subset(novel_I_bed, c(V10 <2))
single_II <- subset(novel_II_bed, c(V10 <2))
single_III <- subset(novel_III_bed, c(V10 <2))
single_intergenic <- subset(intergenic_bed, c(V10 <2))
#Apply a TPM threshold to these transcripts in a tissue specific manner
##first, must attach tissue-specific expression values to these singles that remain
tissue_specific_exp <- read.table("inputs/backmapping_stats/allTissues_isoformTPM", header=T, stringsAsFactors=F)
tissue_specific_intergenic_exp <- read.table("inputs/backmapping_stats/intergenic_allTissues_isoformTPM", header=T, stringsAsFactors=F)

single_known_ncRNA_exp <- merge(single_known,tissue_specific_exp,by.x="transcriptName",by.y="isoformName")
single_known_ncRNA_exp <- single_known_ncRNA_exp[ ,c("transcriptName","BrainStem", "Cerebellum",  "Embryo.ICM", "Embryo.TE",  "Muscle",  "Retina",      "Skin", "SpinalCord")]
single_novel_I_exp <- merge(single_I,tissue_specific_exp,by.x="transcriptName",by.y="isoformName")
single_novel_I_exp <- single_novel_I_exp[ ,c("transcriptName","BrainStem", "Cerebellum",  "Embryo.ICM", "Embryo.TE",  "Muscle",  "Retina",	"Skin",	"SpinalCord")]
single_novel_II_exp <- merge(single_II,tissue_specific_exp,by.x="transcriptName",by.y="isoformName")
single_novel_II_exp <- single_novel_II_exp[ ,c("transcriptName","BrainStem", "Cerebellum",  "Embryo.ICM", "Embryo.TE",	"Muscle",	"Retina",	"Skin",	"SpinalCord")]
single_novel_III_exp <- merge(single_III,tissue_specific_exp,by.x="transcriptName",by.y="isoformName")
single_novel_III_exp <- single_novel_III_exp[ ,c("transcriptName","BrainStem", "Cerebellum",  "Embryo.ICM", "Embryo.TE",	"Muscle",	"Retina",	"Skin",	"SpinalCord")]
single_intergenic_exp <- merge(single_intergenic,tissue_specific_intergenic_exp,by.x="transcriptName",by.y="isoformName")
single_intergenic_exp <- single_intergenic_exp[ ,c("transcriptName","BrainStem", "Cerebellum",  "Embryo.ICM", "Embryo.TE",  "Muscle",	"Retina",	"Skin",	"SpinalCord")]
#now isolate single exons which have at least 5 TPM in any tissue
single_known_ncRNA_2 <- single_known_ncRNA_exp[apply(single_known_ncRNA_exp[-1],1,function(row) {any(row > 5)}),]
single_novel_I_2 <- single_novel_I_exp[apply(single_novel_I_exp[-1],1,function(row) {any(row > 5)}),]
single_novel_II_2 <- single_novel_II_exp[apply(single_novel_II_exp[-1],1,function(row) {any(row > 5)}),]
single_novel_III_2 <- single_novel_III_exp[apply(single_novel_III_exp[-1],1,function(row) {any(row > 5)}),]
single_intergenic_2 <- single_intergenic_exp[apply(single_intergenic_exp[-1],1,function(row) {any(row > 5)}),]
#need to subset the ones removed from this filter
f1_singles_rejects_known <-anti_join(single_known_ncRNA_exp,single_known_ncRNA_2, by="transcriptName")
f1_singles_rejects_I <-anti_join(single_novel_I_exp,single_novel_I_2, by="transcriptName")
f1_singles_rejects_II <-anti_join(single_novel_II_exp,single_novel_II_2, by="transcriptName")
f1_singles_rejects_III <-anti_join(single_novel_III_exp,single_novel_III_2, by="transcriptName")
f1_singles_rejects_intergenic <-anti_join(single_intergenic_exp,single_intergenic_2, by="transcriptName")
#remove the single exon transcripts that did not pass the TPM filter, we move forward with these products
known_ncRNA_f2single <-anti_join(known_ncRNA_f2,f1_singles_rejects_known, by="transcriptName")
novel_I_f2single <-anti_join(novel_I_f2,f1_singles_rejects_I, by="transcriptName")
novel_II_f2single <-anti_join(novel_II_f2,f1_singles_rejects_II, by="transcriptName")
novel_III_f2single <-anti_join(novel_III_f2,f1_singles_rejects_III, by="transcriptName")
novel_intergenic_f2single <-anti_join(intergenic_f2,f1_singles_rejects_intergenic, by="transcriptName")

#convert this output into a bed file
known_ncRNA_f2_bed <- merge(known_ncRNA_f2single, unfiltered_bed, by.x="transcriptName",by.y="V4" )
novel_I_f2_bed <- merge(novel_I_f2single, unfiltered_bed, by.x="transcriptName",by.y="V4" )
novel_II_f2_bed <- merge(novel_II_f2single, unfiltered_bed, by.x="transcriptName",by.y="V4" )
novel_III_f2_bed <- merge(novel_III_f2single, unfiltered_bed, by.x="transcriptName",by.y="V4" )
intergenic_f2_bed <- merge(novel_intergenic_f2single, unfiltered_bed, by.x="transcriptName",by.y="V4" )
#remove non-bed format columns and format properly
known_ncRNA_f2_bed <- known_ncRNA_f2_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")]
names(known_ncRNA_f2_bed)[4]<-paste("V4")
novel_I_f2_bed <- novel_I_f2_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_I_f2_bed)[4]<-paste("V4")
novel_II_f2_bed <- novel_II_f2_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_II_f2_bed)[4]<-paste("V4")
novel_III_f2_bed <- novel_III_f2_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(novel_III_f2_bed)[4]<-paste("V4")
intergenic_f2_bed <- intergenic_f2_bed[ ,c("V1","V2","V3","transcriptName","V5","V6","V7","V8","V9","V10","V11","V12")] 
names(intergenic_f2_bed)[4]<-paste("V4")

#order chr properly
known_ncRNA_f2_bed <- known_ncRNA_f2_bed[with(known_ncRNA_f2_bed, order(V1, V2)), ]
novel_I_f2_bed <- novel_I_f2_bed[with(novel_I_f2_bed, order(V1, V2)), ]
novel_II_f2_bed <- novel_II_f2_bed[with(novel_II_f2_bed, order(V1, V2)), ]
novel_III_f2_bed <- novel_III_f2_bed[with(novel_III_f2_bed, order(V1, V2)), ]
intergenic_f2_bed <- intergenic_f2_bed[with(intergenic_f2_bed, order(V1, V2)), ]
#unfiltered_bed <- unfiltered_bed[with(unfiltered_bed, order(V1, V2)), ]

#write the bedfiles
write.table(known_ncRNA_f2_bed, "known_ncRNA_f2.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_I_f2_bed, "novel_I_f2.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_II_f2_bed, "novel_II_f2.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(novel_III_f2_bed, "novel_III_f2.bed", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(intergenic_f2_bed, "intergenic_f2.bed", row.names=F, col.names=F, quote=F, sep = "\t")
#write.table(unfiltered_bed, "unfiltered.bed", row.names=F, col.names=F, quote=F, sep = "\t")
##use these files to look at intersect with genes (using bedtools -intersect)

#compiling and writing a file for the F1 rejects for each input for single and multiexonic
#f1 for single exonic got switched to f2, hence what looks like mergeing the wrong filters
f1_singles_rejects_known <- merge(f1_singles_rejects_known,overallExpression,by="transcriptName")
f1_singles_known_trunc <- f1_singles_rejects_known[ ,c("transcriptName","length",  "calcTPM")]
F1_known_ncRNA <-rbind(f1_known_rejects,f1_singles_known_trunc)
f1_singles_rejects_I <- merge(f1_singles_rejects_I,overallExpression,by="transcriptName")
f1_singles_I_trunc <- f1_singles_rejects_I[ ,c("transcriptName","length",  "calcTPM")]
F1_novel_I <-rbind(f1_I_rejects,f1_singles_I_trunc)
f1_singles_rejects_II <- merge(f1_singles_rejects_II,overallExpression,by="transcriptName")
f1_singles_II_trunc <- f1_singles_rejects_II[ ,c("transcriptName","length",  "calcTPM")]
F1_novel_II <-rbind(f1_II_rejects,f1_singles_II_trunc)
f1_singles_rejects_III <- merge(f1_singles_rejects_III,overallExpression,by="transcriptName")
f1_singles_III_trunc <- f1_singles_rejects_III[ ,c("transcriptName","length",  "calcTPM")]
F1_novel_III <-rbind(f1_III_rejects,f1_singles_III_trunc)
f1_singles_rejects_intergenic <- merge(f1_singles_rejects_intergenic,overallExpression,by="transcriptName")
f1_singles_intergenic_trunc <- f1_singles_rejects_intergenic[ ,c("transcriptName","length",  "calcTPM")]
F1_intergenic <-rbind(f1_intergenic_rejects,f1_singles_intergenic_trunc)
write.table(F1_known_ncRNA,"F1_known_ncRNA", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(F1_novel_I, "F1_novel_I", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(F1_novel_II, "F1_novel_II", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(F1_novel_III, "F1_novel_III", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(F1_intergenic, "F1_intergenic", row.names=F, col.names=F, quote=F, sep = "\t")
##for F2 rejects
write.table(f2_known_rejects, "F2_known_ncRNA", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(f2_I_rejects, "F2_novel_I", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(f2_II_rejects, "F2_novel_II", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(f2_III_rejects, "F2_novel_III", row.names=F, col.names=F, quote=F, sep = "\t")
write.table(f2_intergenic_rejects, "F2_intergenic", row.names=F, col.names=F, quote=F, sep = "\t")

