##A) cumulative freq vs % seq_identity*%seq_coverage
#genes vs lncRNA against human
setwd("~/lncRNA")
blast_genes <- read.table("coding_hg_trans.outfmt6")
blast_lncRNA <- read.table("lncRNA_hg_trans.outfmt6")
##[3]=%identity,

#make a new column in dataset that is %identity*%coverage
library(plyr)
blast_lncRNA_new<- ddply(blast_lncRNA, c("V1"), summarise,
                         cons=V3*(V4/100))
blast_genes_new<- ddply(blast_genes, c("V1"), summarise,
                        cons=V3*(V4/100))

blast_lncRNA_new <- blast_lncRNA_new[with(blast_lncRNA_new, order(V1, cons, decreasing = T)), ]
blast_lncRNA_new <- blast_lncRNA_new[!duplicated(blast_lncRNA_new$V1),]

blast_genes_new <- blast_genes_new[with(blast_genes_new, order(V1, cons, decreasing = T)), ]
blast_genes_new <- blast_genes_new[!duplicated(blast_genes_new$V1),]

#need to calc or count-freq of each identity*coverage then cumsum(this)=y
blast_lncRNA_freq <- count(blast_lncRNA_new,'cons')
blast_genes_freq <- count(blast_genes_new,'cons')

blast_lncRNA_freq_cumsum <- cbind(blast_lncRNA_freq,cumsum(blast_lncRNA_freq$freq))
blast_genes_freq_cumsum <- cbind(blast_genes_freq,cumsum(blast_genes_freq$freq))

names(blast_lncRNA_freq_cumsum)[3] <-"cumsum" 
names(blast_genes_freq_cumsum)[3] <-"cumsum"
#divide the cumulative sum by the number of rows in the blast_lncRNA and blast_genes
#to get cumulative relative freq
#blast_lncRNA_freq_cumsum["cum_rel"]<-(blast_lncRNA_freq_cumsum$cumsum/dim(blast_lncRNA_new)[1])
#blast_genes_freq_cumsum["cum_rel"]<-(blast_genes_freq_cumsum$cumsum/dim(blast_genes_new)[1])

#blast_lncRNA_gene_cumsum_rel <- rbind(data.frame(id="lncRNA",blast_lncRNA_freq_cumsum),
#                                      data.frame(id="PCG",blast_genes_freq_cumsum))
#write.table(blast_lncRNA_gene_cumsum_rel, "blast_lncRNA_gene_cumsum_rel.txt")
# plot
#require(ggplot2)
#pdf("Fig3A.pdf")
#ggplot(blast_lncRNA_gene_cumsum_rel, aes(x=cons, y=cum_rel, colour=id, group=id)) + 
#  geom_line() + xlab("Blast %identity* Blast % coverage") + 
#  ylab("cumulative frequency") + scale_color_discrete(name="Type of annotation")
#dev.off()

## Re-run using total transcript count (replace line 32 - 44)
nc_count=length(readLines("lncRNA_final.bed"))
cod_count=length(readLines("refined_codingRNA.bed"))
blast_lncRNA_freq_cumsum[1,2]= blast_lncRNA_freq_cumsum[1,2] + (nc_count - dim(blast_lncRNA_new)[1])
blast_lncRNA_freq_cumsum[,3]= blast_lncRNA_freq_cumsum[,3] + (nc_count - dim(blast_lncRNA_new)[1])
blast_genes_freq_cumsum[1,2]= blast_genes_freq_cumsum[1,2] + (cod_count - dim(blast_genes_new)[1]) 
blast_genes_freq_cumsum[,3]= blast_genes_freq_cumsum[,3] + (cod_count - dim(blast_genes_new)[1])

blast_lncRNA_freq_cumsum["cum_rel"]<-(blast_lncRNA_freq_cumsum$cumsum/nc_count) 
blast_genes_freq_cumsum["cum_rel"]<-(blast_genes_freq_cumsum$cumsum/cod_count) 

blast_lncRNA_gene_cumsum_rel <- rbind(data.frame(id="lncRNA",blast_lncRNA_freq_cumsum),
                                      data.frame(id="PCG",blast_genes_freq_cumsum))
write.table(blast_lncRNA_gene_cumsum_rel, "blast_lncRNA_gene_cumsum_rel_v2.txt")
# plot
require(ggplot2)
pdf("Fig3A_v2.pdf")
ggplot(blast_lncRNA_gene_cumsum_rel, aes(x=cons, y=cum_rel, colour=id, group=id)) +
  geom_line() + xlab("Blast %identity* Blast % coverage") +
  ylab("cumulative frequency") + scale_color_discrete(name="Type of annotation")
dev.off()

##for promoters and down of lncRNA and genes
blast_genes_promoter <- read.table("genepromoter_human_nt.outfmt6")
blast_genes_down <- read.table("genedown_human_nt.outfmt6")
blast_lncRNA_promoter <- read.table("lncRNApromoter_human_nt.outfmt6")
blast_lncRNA_down <- read.table("lncRNAdown_human_nt.outfmt6")
blast_lncRNAg <- read.table("lncRNA_hg_genome.outfmt6")

library(plyr)
blast_lncRNA_promoter_new <- ddply(blast_lncRNA_promoter, c("V1"), summarise, cons=V3*(V4/100))
blast_lncRNA_down_new <- ddply(blast_lncRNA_down, c("V1"), summarise, cons=V3*(V4/100))
blast_genes_promoter_new <- ddply(blast_genes_promoter, c("V1"), summarise, cons=V3*(V4/100))
blast_genes_down_new <- ddply(blast_genes_down, c("V1"), summarise, cons=V3*(V4/100))
blast_lncRNAg_new <- ddply(blast_lncRNAg, c("V1"), summarise, cons=V3*(V4/100))

blast_lncRNA_promoter_new <- blast_lncRNA_promoter_new[with(blast_lncRNA_promoter_new, order(V1, cons, decreasing = T)), ]
blast_lncRNA_promoter_new <- blast_lncRNA_promoter_new[!duplicated(blast_lncRNA_promoter_new$V1),]

blast_lncRNA_down_new <- blast_lncRNA_down_new[with(blast_lncRNA_down_new, order(V1, cons, decreasing = T)), ]
blast_lncRNA_down_new <- blast_lncRNA_down_new[!duplicated(blast_lncRNA_down_new$V1),]

blast_genes_promoter_new <- blast_genes_promoter_new[with(blast_genes_promoter_new, order(V1, cons, decreasing = T)), ]
blast_genes_promoter_new <- blast_genes_promoter_new[!duplicated(blast_genes_promoter_new$V1),]

blast_genes_down_new <- blast_genes_down_new[with(blast_genes_down_new, order(V1, cons, decreasing = T)), ]
blast_genes_down_new <- blast_genes_down_new[!duplicated(blast_genes_down_new$V1),]

blast_lncRNAg_new <- blast_lncRNAg_new[with(blast_lncRNAg_new, order(V1, cons, decreasing = T)), ]
blast_lncRNAg_new <- blast_lncRNAg_new[!duplicated(blast_lncRNAg_new$V1),]

blast_lncRNA_promoter_freq <- count(blast_lncRNA_promoter_new,'cons')
blast_lncRNA_down_freq <- count(blast_lncRNA_down_new,'cons')
blast_genes_promoter_freq <- count(blast_genes_promoter_new,'cons')
blast_genes_down_freq <- count(blast_genes_down_new,'cons')
blast_lncRNAg_freq <- count(blast_lncRNAg_new,'cons')

#blank_count=data.frame(cons=0,freq=0)
#blast_lncRNA_promoter_freq=rbind(blank_count,blast_lncRNA_promoter_freq)
#blast_lncRNA_down_freq=rbind(blank_count,blast_lncRNA_down_freq)
#blast_genes_promoter_freq=rbind(blank_count,blast_genes_promoter_freq)
#blast_genes_down_freq=rbind(blank_count,blast_genes_down_freq)

blast_lncRNA_promoter_cumsum <- cbind(blast_lncRNA_promoter_freq,cumsum(blast_lncRNA_promoter_freq$freq))
blast_lncRNA_down_cumsum <- cbind(blast_lncRNA_down_freq,cumsum(blast_lncRNA_down_freq$freq))
blast_genes_promoter_cumsum <- cbind(blast_genes_promoter_freq,cumsum(blast_genes_promoter_freq$freq))
blast_genes_down_cumsum <- cbind(blast_genes_down_freq,cumsum(blast_genes_down_freq$freq))
blast_lncRNAg_cumsum <- cbind(blast_lncRNAg_freq,cumsum(blast_lncRNAg_freq$freq))

names(blast_lncRNA_promoter_cumsum)[3] <-"cumsum" 
names(blast_lncRNA_down_cumsum)[3] <-"cumsum" 
names(blast_genes_promoter_cumsum)[3] <-"cumsum" 
names(blast_genes_down_cumsum)[3] <-"cumsum" 
names(blast_lncRNAg_cumsum)[3] <-"cumsum"

#divide the cumulative sum by the number of rows in the blast_lncRNA and blast_genes
#to get cumulative relative freq
#blast_lncRNA_promoter_cumsum["cum_rel"]<-(blast_lncRNA_promoter_cumsum$cumsum/dim(blast_lncRNA_promoter_new)[1])
#blast_lncRNA_down_cumsum["cum_rel"]<-(blast_lncRNA_down_cumsum$cumsum/dim(blast_lncRNA_down_new)[1])
#blast_genes_promoter_cumsum["cum_rel"]<-(blast_genes_promoter_cumsum$cumsum/dim(blast_genes_promoter_new)[1])
#blast_genes_down_cumsum["cum_rel"]<-(blast_genes_down_cumsum$cumsum/dim(blast_genes_down_new)[1])

#blast_lncRNA_gene_PD_cumsum_rel <- rbind(data.frame(id="lncRNA_promoters",blast_lncRNA_promoter_cumsum),
#                                         data.frame(id="lncRNA_down",blast_lncRNA_down_cumsum),
#                                         data.frame(id="gene_promoters",blast_genes_promoter_cumsum),
#                                         data.frame(id="gene_down",blast_genes_down_cumsum))
#write.table(blast_lncRNA_gene_PD_cumsum_rel,"blast_lncRNA_gene_PD_cumsum_rel.txt")                                     
# plot
#require(ggplot2)
#pdf("Fig3B.pdf")
#ggplot(blast_lncRNA_gene_PD_cumsum_rel, aes(x=cons, y=cum_rel, colour=id, group=id)) + 
#  geom_line(position=position_dodge(width=0.2)) + xlab("Blast %identity* Blast % coverage") + 
#  geom_line() + xlab("Blast %identity* Blast % coverage") +
#  ylab("cumulative frequency") + scale_color_discrete(name="Type of annotation")
#dev.off()

## Re-run using total transcript count (replace line 112 - 129)
ncPro_count=length(readLines("lncRNA_promoters_uniq.bed"))
ncDown_count=length(readLines("lncRNA_down_uniq.bed"))
codPro_count=length(readLines("refined_codingRNA_5_uniq.bed"))
codDown_count=length(readLines("refined_codingRNA_3_uniq.bed"))
nc_count=length(readLines("lncRNA_final.bed"))

blast_lncRNA_promoter_cumsum[1,2]= blast_lncRNA_promoter_cumsum[1,2] + (ncPro_count - dim(blast_lncRNA_promoter_new)[1])
blast_lncRNA_promoter_cumsum[,3]= blast_lncRNA_promoter_cumsum[,3] + (ncPro_count - dim(blast_lncRNA_promoter_new)[1])
blast_lncRNA_down_cumsum[1,2]= blast_lncRNA_down_cumsum[1,2] + (ncDown_count - dim(blast_lncRNA_down_new)[1])
blast_lncRNA_down_cumsum[,3]= blast_lncRNA_down_cumsum[,3] + (ncDown_count - dim(blast_lncRNA_down_new)[1])
blast_genes_promoter_cumsum[1,2]= blast_genes_promoter_cumsum[1,2] + (codPro_count - dim(blast_genes_promoter_new)[1])
blast_genes_promoter_cumsum[,3]= blast_genes_promoter_cumsum[,3] + (codPro_count - dim(blast_genes_promoter_new)[1])
blast_genes_down_cumsum[1,2]= blast_genes_down_cumsum[1,2] + (codDown_count - dim(blast_genes_down_new)[1]) 
blast_genes_down_cumsum[,3]= blast_genes_down_cumsum[,3] + (codDown_count - dim(blast_genes_down_new)[1])
blast_lncRNAg_cumsum[1,2]= blast_lncRNAg_cumsum[1,2] + (nc_count - dim(blast_lncRNAg_new)[1])
blast_lncRNAg_cumsum[,3]= blast_lncRNAg_cumsum[,3] + (nc_count - dim(blast_lncRNAg_new)[1])

blast_lncRNA_promoter_cumsum["cum_rel"]<-(blast_lncRNA_promoter_cumsum$cumsum/ncPro_count)
blast_lncRNA_down_cumsum["cum_rel"]<-(blast_lncRNA_down_cumsum$cumsum/ncDown_count)
blast_genes_promoter_cumsum["cum_rel"]<-(blast_genes_promoter_cumsum$cumsum/codPro_count)
blast_genes_down_cumsum["cum_rel"]<-(blast_genes_down_cumsum$cumsum/codDown_count)
blast_lncRNAg_cumsum["cum_rel"]<-(blast_lncRNAg_cumsum$cumsum/nc_count)

blast_lncRNA_gene_PD_cumsum_rel <- rbind(data.frame(id="lncRNA_promoters",blast_lncRNA_promoter_cumsum),
                                         data.frame(id="lncRNA_down",blast_lncRNA_down_cumsum),
                                         data.frame(id="gene_promoters",blast_genes_promoter_cumsum),
                                         data.frame(id="gene_down",blast_genes_down_cumsum))

write.table(blast_lncRNA_gene_PD_cumsum_rel,"blast_lncRNA_gene_PD_cumsum_rel_v2.txt")
# plot
require(ggplot2)
pdf("Fig3B_v2.pdf")
ggplot(blast_lncRNA_gene_PD_cumsum_rel, aes(x=cons, y=cum_rel, colour=id, group=id)) +
#  geom_line(position=position_dodge(width=0.2)) + xlab("Blast %identity* Blast % coverage") + 
  geom_line() + xlab("Blast %identity* Blast % coverage") +
  ylab("cumulative frequency") + scale_color_discrete(name="Type of annotation")
dev.off()

blast_lncRNA_gene_PD_cumsum_rel <- rbind(data.frame(id="lncRNA",blast_lncRNAg_cumsum),
                                         data.frame(id="lncRNA_promoters",blast_lncRNA_promoter_cumsum),
                                         data.frame(id="lncRNA_down",blast_lncRNA_down_cumsum),
                                         data.frame(id="gene_promoters",blast_genes_promoter_cumsum),
                                         data.frame(id="gene_down",blast_genes_down_cumsum))

write.table(blast_lncRNA_gene_PD_cumsum_rel,"blast_lncRNA_gene_PD_cumsum_rel_v3.txt")
# plot
require(ggplot2)
pdf("Fig3B_v3.pdf")
ggplot(blast_lncRNA_gene_PD_cumsum_rel, aes(x=cons, y=cum_rel, colour=id, group=id)) +
#  geom_line(position=position_dodge(width=0.2)) + xlab("Blast %identity* Blast % coverage") + 
  geom_line() + xlab("Blast %identity* Blast % coverage") +
  ylab("cumulative frequency") + scale_color_discrete(name="Type of annotation")
dev.off()

