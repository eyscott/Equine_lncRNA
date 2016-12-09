setwd("~/Desktop/lncRNA")
novel_I_GC <- read.table("novel_I_GC.bed", header=F, stringsAsFactors=F)
novel_II_GC <- read.table("novel_II_GC.bed", header=F, stringsAsFactors=F)
novel_III_GC <- read.table("novel_III_GC.bed", header=F, stringsAsFactors=F)
intergenic_GC <- read.table("intergenic_GC.bed", header=F, stringsAsFactors=F)
known_GC <- read.table("known_GC.bed", header=F, stringsAsFactors=F)

#calculating mean GC%
meanGC_novel_I<-mean(novel_I_GC[["V14"]])
meanGC_novel_II<-mean(novel_II_GC[["V14"]])
meanGC_novel_III<-mean(novel_III_GC[["V14"]])
meanGC_intergenic<-mean(intergenic_GC[["V14"]])
meanGC_known<-mean(known_GC[["V14"]])
