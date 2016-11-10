#lncRNA quality stats from Transdecoder
setwd("~/Desktop/lncRNA")
#interpretting the Transdecoder stats
#getting sequence quality
setwd("~/Desktop/lncRNA/intergenic.fa.transdecoder_dir")
intergenic_seqQ <- read.table("hexamer.scores", header=F, stringsAsFactors=F)
setwd("~/Desktop/lncRNA/novel_I.fa.transdecoder_dir")
novel_I_seqQ <- read.table("hexamer.scores", header=F, stringsAsFactors=F)
setwd("~/Desktop/lncRNA/novel_II.fa.transdecoder_dir")
novel_II_seqQ <- read.table("hexamer.scores", header=F, stringsAsFactors=F)
setwd("~/Desktop/lncRNA/novel_III.fa.transdecoder_dir")
novel_III_seqQ <- read.table("hexamer.scores", header=F, stringsAsFactors=F)
#Getting ORF quality
setwd("~/Desktop/lncRNA/intergenic.fa.transdecoder_dir")
intergenic_ORF <- read.table("longest_orfs.cds.scores", header=F, stringsAsFactors=F)
setwd("~/Desktop/lncRNA/novel_I.fa.transdecoder_dir")
novel_I_ORF <- read.table("longest_orfs.cds.scores", header=F, stringsAsFactors=F)
setwd("~/Desktop/lncRNA/novel_II.fa.transdecoder_dir")
novel_II_ORF <- read.table("longest_orfs.cds.scores", header=F, stringsAsFactors=F)
setwd("~/Desktop/lncRNA/novel_III.fa.transdecoder_dir")
novel_III_ORF <- read.table("longest_orfs.cds.scores", header=F, stringsAsFactors=F)

#calculating seq quality based on hexamer, k-mer
meanQ_novel_I<-mean(novel_I_seqQ[["V2"]])
meanQ_novel_II<-mean(novel_II_seqQ[["V2"]])
meanQ_novel_III<-mean(novel_III_seqQ[["V2"]])
meanQ_intergenic<-mean(intergenic_seqQ[["V2"]])

#calculating ORF quality based on likely hood ratio test
meanORF_novel_I<-mean(novel_I_ORF[["V2"]])
meanORF_novel_II<-mean(novel_II_ORF[["V2"]])
meanORF_novel_III<-mean(novel_III_ORF[["V2"]])
meanORF_intergenic<-mean(intergenic_ORF[["V2"]])

#compiling GC% after using bedtools nuc (V13=AT%, V14=GC%)
setwd("~/Desktop/lncRNA")
novel_I_GC <- read.table("novel_I_GC.bed", header=F, stringsAsFactors=F)
novel_II_GC <- read.table("novel_II_GC.bed", header=F, stringsAsFactors=F)
novel_III_GC <- read.table("novel_III_GC.bed", header=F, stringsAsFactors=F)
intergenic_GC <- read.table("intergenic_GC.bed", header=F, stringsAsFactors=F)
#calculating mean GC%
meanGC_novel_I<-mean(novel_I_GC[["V14"]])
meanGC_novel_II<-mean(novel_II_GC[["V14"]])
meanGC_novel_III<-mean(novel_III_GC[["V14"]])
meanGC_intergenic<-mean(intergenic_GC[["V14"]])

