#Figure 3: tissue-specific expression of lncRNA
#A) lncRNA vs annotated genes vs uncertainty
#Use the mergedTrans.bed without the novel I,II,III
setwd("~/Dropbox/Horse_Transcriptome/downloads")
mergedTrans <- read.table("allTissues_BED/mergedTrans.bed", header=F, stringsAsFactors=F)
mergedTrans <- mergedTrans[ ,c("V1","V2","V3","V4")]
setwd("~/Desktop/lncRNA")
all <- read.table("lncRNA_final21032", header=F, stringsAsFactors=F)
#make sure the mergedtrans has none of the novel I,II,III
novels <- read.table("all_cats_PandnoP", header=F, stringsAsFactors=F)
novels <- novels[ ,c("V2","V3","V4","V5")]
names(novels)<-names(mergedTrans)
library(dplyr)
mergedTrans_noNovel <- anti_join(mergedTrans,novels)
#join with expression data
setwd("~/Dropbox/Horse_Transcriptome/downloads")
tissue_specific_intergenic_exp <- read.table("intergenic_trans/allTissues_isoformTPM", header=T, stringsAsFactors=F)
tissue_specific_exp <- read.table("backmapping_stats/allTissues_isoformTPM", header=T, stringsAsFactors=F)
#remove mt entries
tissue_specific_exp <- tissue_specific_exp[-c(1,2),]
rownames(tissue_specific_exp) <- c()
#combine
tissue_specific_all_exp <- rbind(tissue_specific_intergenic_exp,tissue_specific_exp)
##merge with 
setwd("~/Desktop/lncRNA")
not_annotated <- read.table("nonannotated", header=T, stringsAsFactors=F)
not_annotated <- merge(not_annotated,tissue_specific_all_exp, by.x="transcriptName",by.y="isoformName")
lncRNA <- subset(not_annotated, V1 %in% c("novel_I_lncRNA","novel_II_lncRNA","novel_II_lncRNA","intergenic_lncRNA"))
uncertain <- anti_join(not_annotated,lncRNA, by.x="transcripName",by.y="transcripName")
lncRNA <- lncRNA[ ,c("transcriptName","BrainStem", "Cerebellum",  "Embryo.ICM", "Embryo.TE",  "Muscle",  "Retina",	"Skin",	"SpinalCord")]
rownames(lncRNA) <- c()
uncertain <- uncertain[ ,c("transcriptName","BrainStem", "Cerebellum",  "Embryo.ICM", "Embryo.TE",  "Muscle",  "Retina",  "Skin",	"SpinalCord")]
#calculate cumulative TPM for pie graph
cumulative_uncertain_TPM_BS<-sum(uncertain[["BrainStem"]])
cumulative_uncertain_TPM_C<-sum(uncertain[["Cerebellum"]])
cumulative_uncertain_TPM_EICM<-sum(uncertain[["Embryo.ICM"]])
cumulative_uncertain_TPM_ETE<-sum(uncertain[["Embryo.TE"]])
cumulative_uncertain_TPM_M<-sum(uncertain[["Muscle"]])
cumulative_uncertain_TPM_R<-sum(uncertain[["Retina"]])
cumulative_uncertain_TPM_S<-sum(uncertain[["Skin"]])
cumulative_uncertain_TPM_SC<-sum(uncertain[["SpinalCord"]])

mergedTrans_noNovel_exp <- merge(mergedTrans_noNovel,tissue_specific_all_exp, by.x="V4",by.y="isoformName")
mergedTrans_noNovel_exp <- mergedTrans_noNovel_exp[ , c("V4","BrainStem", "Cerebellum",  "Embryo.ICM", "Embryo.TE",  "Muscle",	"Retina",	"Skin",	"SpinalCord")] 
names(mergedTrans_noNovel_exp)[1]<-paste("transcriptName")
#calculate cumulative TPM for pie graph
cumulative_annotated_TPM_BS<-sum(mergedTrans_noNovel_exp[["BrainStem"]])
cumulative_annotated_TPM_C<-sum(mergedTrans_noNovel_exp[["Cerebellum"]])
cumulative_annotated_TPM_EICM<-sum(mergedTrans_noNovel_exp[["Embryo.ICM"]])
cumulative_annotated_TPM_ETE<-sum(mergedTrans_noNovel_exp[["Embryo.TE"]])
cumulative_annotated_TPM_M<-sum(mergedTrans_noNovel_exp[["Muscle"]])
cumulative_annotated_TPM_R<-sum(mergedTrans_noNovel_exp[["Retina"]])
cumulative_annotated_TPM_S<-sum(mergedTrans_noNovel_exp[["Skin"]])
cumulative_annotated_TPM_SC<-sum(mergedTrans_noNovel_exp[["SpinalCord"]])

#combine all
total_exp <- rbind(data.frame(id="lncRNA",lncRNA),
                   data.frame(id="uncertain",uncertain),
                   data.frame(id="annotated",mergedTrans_noNovel_exp))
BrainStem <- total_exp[ ,c("id","BrainStem")]
Cerebellum <- total_exp[ ,c("id","Cerebellum")]
Embryo.ICM <- total_exp[ ,c("id","Embryo.ICM")]
Embryo.TE <- total_exp[ ,c("id","Embryo.TE")]
Muscle <- total_exp[ ,c("id","Muscle")]
Retina <- total_exp[ ,c("id","Retina")]
Skin <- total_exp[ ,c("id","Skin")]
SpinalCord <- total_exp[ ,c("id","SpinalCord")]
# save each of these for distribution plot
write.table(total_exp, "total_tissue_all_TPM", row.names=F, col.names=T, sep = "\t")
#count number of lncRNA/tissue
BS_lncRNA<-subset(BrainStem,id %in% "lncRNA")#18025
C_lncRNA<-subset(Cerebellum,id %in% "lncRNA")#18025
EICM_lncRNA<-subset(Embryo.ICM,id %in% "lncRNA")
ETE_lncRNA<-subset(Embryo.TE,id %in% "lncRNA")
M_lncRNA<-subset(Muscle,id %in% "lncRNA")
R_lncRNA<-subset(Retina,id %in% "lncRNA")
S_lncRNA<-subset(Skin,id %in% "lncRNA")
SC_lncRNA<-subset(SpinalCord,id %in% "lncRNA")
#make a cutoff of 0.1 TPM
BS_lncRNA_cut<-subset(BS_lncRNA,BrainStem > 0.1)#13308,genes:34423
rownames(BS_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_BS<-sum(BS_lncRNA_cut[["BrainStem"]])
C_lncRNA_cut<-subset(C_lncRNA,Cerebellum > 0.1)#14878,genes:35784
rownames(C_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_C<-sum(C_lncRNA_cut[["Cerebellum"]])
EICM_lncRNA_cut<-subset(EICM_lncRNA,Embryo.ICM > 0.1)#15640,genes:33566
rownames(EICM_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_EICM<-sum(EICM_lncRNA_cut[["Embryo.ICM"]])
ETE_lncRNA_cut<-subset(ETE_lncRNA,Embryo.TE > 0.1)#11209,genes:31633
rownames(ETE_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_ETE<-sum(ETE_lncRNA_cut[["Embryo.TE"]])
M_lncRNA_cut<-subset(M_lncRNA,Muscle > 0.1)#2848,,genes:29203
rownames(M_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_M<-sum(M_lncRNA_cut[["Muscle"]])
R_lncRNA_cut<-subset(R_lncRNA,Retina > 0.1)#3934,genes:26419
rownames(R_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_R<-sum(R_lncRNA_cut[["Retina"]])
S_lncRNA_cut<-subset(S_lncRNA,Skin > 0.1)#5417,genes:29674
rownames(S_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_S<-sum(S_lncRNA_cut[["Skin"]])
SC_lncRNA_cut<-subset(SC_lncRNA,SpinalCord > 0.1)#12736,genes:34636
rownames(SC_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_SC<-sum(SC_lncRNA_cut[["SpinalCord"]])
#making lists for pie
BrainStem = list(nv(c(cumulative_lncRNA_TPM_BS,cumulative_uncertain_TPM_BS,cumulative_annotated_TPM_BS),c('lncRNA','uncertain','annotated')))
Cerebellum = list(nv(c(cumulative_lncRNA_TPM_C,cumulative_uncertain_TPM_C,cumulative_annotated_TPM_C),c('lncRNA','uncertain','annotated')))
Embryo.ICM = list(nv(c(cumulative_lncRNA_TPM_EICM,cumulative_uncertain_TPM_EICM,cumulative_annotated_TPM_EICM),c('lncRNA','uncertain','annotated')))
Embryo.TE = list(nv(c(cumulative_lncRNA_TPM_ETE,cumulative_uncertain_TPM_ETE,cumulative_annotated_TPM_ETE),c('lncRNA','uncertain','annotated')))
Muscle = list(nv(c(cumulative_lncRNA_TPM_M,cumulative_uncertain_TPM_M,cumulative_annotated_TPM_M),c('lncRNA','uncertain','annotated')))
Retina = list(nv(c(cumulative_lncRNA_TPM_R,cumulative_uncertain_TPM_R,cumulative_annotated_TPM_R),c('lncRNA','uncertain','annotated')))
Skin = list(nv(c(cumulative_lncRNA_TPM_S,cumulative_uncertain_TPM_S,cumulative_annotated_TPM_S),c('lncRNA','uncertain','annotated')))
SpinalCord = list(nv(c(cumulative_lncRNA_TPM_SC,cumulative_uncertain_TPM_SC,cumulative_annotated_TPM_SC),c('lncRNA','uncertain','annotated')))

BrainStem = list(c(cumulative_lncRNA_TPM_BS,cumulative_uncertain_TPM_BS,cumulative_annotated_TPM_BS))
Cerebellum = list(c(cumulative_lncRNA_TPM_C,cumulative_uncertain_TPM_C,cumulative_annotated_TPM_C))
Embryo.ICM = list(c(cumulative_lncRNA_TPM_EICM,cumulative_uncertain_TPM_EICM,cumulative_annotated_TPM_EICM))
Embryo.TE = list(c(cumulative_lncRNA_TPM_ETE,cumulative_uncertain_TPM_ETE,cumulative_annotated_TPM_ETE))
Muscle = list(c(cumulative_lncRNA_TPM_M,cumulative_uncertain_TPM_M,cumulative_annotated_TPM_M))
Retina = list(c(cumulative_lncRNA_TPM_R,cumulative_uncertain_TPM_R,cumulative_annotated_TPM_R))
Skin = list(c(cumulative_lncRNA_TPM_S,cumulative_uncertain_TPM_S,cumulative_annotated_TPM_S))
SpinalCord = list(c(cumulative_lncRNA_TPM_SC,cumulative_uncertain_TPM_SC,cumulative_annotated_TPM_SC))

pie_input <- list(
  BrainStem,
  Cerebellum,
  Embryo.ICM,
  Embryo.TE,
  Muscle,
  Retina,
  Skin,
  SpinalCord)

pie.list <- lapply(pie_input, table)
#got gene numbers tissueSpecificSummary
##scatter plot with pie dots
library(caroline)
#calc cumulative sum of lncRNA,uncertain and annotated then input manually?
pies(
  list(
    BrainStem=nv(c(92915.81,946066.50,804372.28),c('lncRNA','uncertain','annotated')),
    Cerebellum=nv(c(145327.4,1163073.2,919594.8 ),c('lncRNA','uncertain','annotated')),
    Embryo.ICM=nv(c(70288.96,1423981.62,618838.78),c('lncRNA','uncertain','annotated')),
    Embryo.TE=nv(c(47578.09,1543808.55,590603.21 ),c('lncRNA','uncertain','annotated')),
    Muscle=nv(c(6726.063,1072639.073,803600.748),c('lncRNA','uncertain','annotated')),
    Retina=nv(c(92837.49,1256694.18,900396.30),c('lncRNA','uncertain','annotated')),
    Skin=nv(c(49757.1,952356.4,899838.1 ),c('lncRNA','uncertain','annotated')),
    SpinalCord=nv(c(130291.6,973360.8 ,860672.0 ),c('lncRNA','uncertain','annotated'))),
  x0=c(34423,35784,33566,31633,29203,26419,29674,34636),
  y0=c(13308,14878,15640,11209,2848,3934,5417,12736),
  radii=6, border=c('red','red','black','black','yellow','yellow','yellow','red'),edges = 8000)
#Pies only including genes and lncRNA
pies(
  list(
    BrainStem=nv(c(92915.81,804372.28),c('lncRNA','annotated')),
    Cerebellum=nv(c(145327.4,919594.8 ),c('lncRNA','annotated')),
    Embryo.ICM=nv(c(70288.96,618838.78),c('lncRNA','annotated')),
    Embryo.TE=nv(c(47578.09,590603.21 ),c('lncRNA','annotated')),
    Muscle=nv(c(6726.063,803600.748),c('lncRNA','annotated')),
    Retina=nv(c(92837.49,900396.30),c('lncRNA','annotated')),
    Skin=nv(c(49757.1,899838.1 ),c('lncRNA','annotated')),
    SpinalCord=nv(c(130291.6,860672.0 ),c('lncRNA','annotated'))),
  x0=c(34423,35784,33566,31633,29203,26419,29674,34636),
  y0=c(13308,14878,15640,11209,2848,3934,5417,12736),
  radii=6, border=c('red','red','black','black','yellow','yellow','yellow','red'),edges = 8000)



#B) tissue-specific heatmap
all_exp <- read.table("lncRNA_tissue_exp", header=T, stringsAsFactors=F)
#Melt data for manipulation
melted <- melt(all_exp, id.vars=c("V4"))
#Calculate the sum(TPM) and STDEV of each gene per tissue
library(plyr)
melted_new<- ddply(melted, c("V4"), summarise,
                   sum = sum(value), sd = sd(value),
                   sem = sd(value)/sqrt(length(value)))
#Add the column of sum and sd to the original TPM values table
complete <- merge(all_exp,melted_new,by="V4")
#Subset data based on if sum>50 and sd>50
sub <- subset(complete, c(sum > 20 & sd > 10))
#make row.names the geneName
rownames(sub)<-sub$V4
rownames(all_exp)<-all_exp$V4
# making the matrix for the heatmap
#disable scientific notation so no "e+/-"
options("scipen"=100, "digits"=4)
datanumbers <- data.matrix(all_exp[,2:9])
datanumbers_smalls <- data.matrix(sub[,2:9])
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("Blue", "white", "Red"))(n = 18)
#making the heatmap
library(gplots)
library(RColorBrewer)
library(svDialogs)

#Another way to manually cluster to your liking...allows use you to choose which
#correlations work best with your data depending on linear or monotonic relationship
#of variables
# Row clustering...pearson seems to work best here 
hr <- hclust(as.dist(1-cor(t(datanumbers_smalls), method="pearson")),
             method="average")
# Column clustering...spearman seems to work best here 
hc <- hclust(as.dist(1-cor(datanumbers_smalls, method="spearman")), method="average")

## Plot heatmap
setwd("~/Desktop/lncRNA")
par(mar=c(7,4,4,2)+0.1) 
png(filename='heatmap_manual_SEM_lncRNA.png', width=800, height=750)
col_breaks <- c(1:10,20,30,40,50,60,70,80,90,100)
heatmap.2(datanumbers_smalls,    # data matrix
          #cellnote = mat_data,  # same data set for cell labels
          #main = "Rank", # heat map title
          #notecex=0.4,
          #cexRow=0.6,
          cexCol=2,
          scale="none",
          #notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,12),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="column",     # only draw a row dendrogram
          Colv=as.dendrogram(hc),
          Rowv=as.dendrogram(hr),
          hclustfun = hclust,
          labRow = NULL,
          key.xlab = "TPM",
          key.title = NULL)            # turn off column clustering
graphics.off()  # close the PNG device
## Return matrix with row/column sorting as in heatmap
write.csv(datanumbers[rev(hr$labels[hr$order]), hc$labels[hc$order]],"hmap_order.csv")
#match the XLOC_* names with real gene names and add other annotation information
#So you can see gene name and chr location,strand, and how it compares in other databases
hmap_order <- read.csv("hmap_order.csv", header=T)
annotated <- read.table("RNAseqSupTrans.merge.reduced", header=T)
m <- merge(hmap_order, annotated, by.x="X",by.y="gene.ID", sort=F)
write.csv(m,"hmap_order_names.csv")

#C) unique vs absent
setwd("~/Dropbox/lncRNA/uniq_exp")
library(ggplot2)
library(reshape2)
library(dplyr)
library(plyr)

###plotting figures with varying threshold for absent vs unique lncRNA
data_0.1<-read.table("tissueSpecificSummary_cutoff.0.1")
data_changed_0.1 <- cbind(as.data.frame(data_0.1[1:4,]),as.data.frame(data_0.1[5:8,]),as.data.frame(data_0.1[9:12,]),as.data.frame(data_0.1[13:16,]),as.data.frame(data_0.1[17:20,]),as.data.frame(data_0.1[21:24,]),as.data.frame(data_0.1[25:28,]),as.data.frame(data_0.1[29:32,]))   
data_changed_0.1 <- sapply(data_changed_0.1, as.character)
colnames(data_changed_0.1) <- data_changed_0.1[1,]
data_changed_0.1 <- as.data.frame(data_changed_0.1[-1,])
rownames(data_changed_0.1) <- c("total","unique_lncRNA","not_unique_lncRNA")
data_0.1 <-as.data.frame(t(data_changed_0.1))
write.table(data_0.1, "tissue_lncRNA_0.1.txt")

data_0.1 <- read.table("tissue_lncRNA_0.1.txt",stringsAsFactors=FALSE)
data_0.1$not_unique_lncRNA <- data_0.1$not_unique_lncRNA * -1


Absent_U_data_0.1 <- as.data.frame(data_0.1$not_unique_lncRNA)
rownames(Absent_U_data_0.1) <- rownames(data_0.1) 
U_data_0.1 <- as.data.frame(data_0.1$unique_lncRNA)
rownames(U_data_0.1) <- rownames(data_0.1)
# to get cumulative TPM of those not expressed
BrainStem <- read.table("TPM/BrainStem.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
Cerebellum <- read.table("TPM/Cerebellum.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
Embryo.ICM <- read.table("TPM/Embryo.ICM.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
Embryo.TE <- read.table("TPM/Embryo.TE.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
Muscle <- read.table("TPM/Muscle.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
Retina <- read.table("TPM/Retina.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
Skin <- read.table("TPM/Skin.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
SpinalCord <- read.table("TPM/SpinalCord.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
all_unique_exp <- rbind(BrainStem,Cerebellum,Embryo.ICM,Embryo.TE,Muscle,Retina,Skin,SpinalCord)
names(all_unique_exp)<-c("id","BrainStem","Cerebellum","Embryo.ICM","Embryo.TE","Muscle","Retina","Skin","SpinalCord")
melt_all_unique_exp <- melt(all_unique_exp,id.vars=c("id"))
melt_all_unique_exp_stats<- ddply(melt_all_unique_exp, c("variable"), summarise,
                      sum = sum(value))

ggplot() +
  geom_bar(data=U_data_0.1, aes(x=rownames(data_0.1),y=data_0.1$unique_lncRNA,color="aliceblue"), stat="identity") +
  geom_bar(data=Absent_U_data_0.1, aes(x=rownames(data_0.1),y=data_0.1$not_unique_lncRNA,color="red"),stat = "identity") + 
  ylab("Number of lncRNA") + scale_color_discrete(name="Unique lncRNA",
                                                    labels=c("present","absent")) + xlab("Tissue") +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(colour="black", size = 9)) +
  theme(axis.title = element_text(colour="black", size = 14)) +
  geom_line(data=melt_all_unique_exp_stats, aes(x=variable,y=sum / 5, group=1),colour="green")

