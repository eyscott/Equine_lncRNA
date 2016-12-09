## Figure 4: tissue-specific expression of lncRNA
#A) lncRNA vs annotated genes
#Use the mergedTrans.bed without the novel I,II,III
setwd("~/lncRNA")
refined_nolncRNA <- read.table("refined_nolncRNA.bed", header=F, stringsAsFactors=F)
all <- read.table("lncRNA_final.bed", header=F, stringsAsFactors=F)

library(plyr)
library(dplyr)
#join with expression data
tissue_specific_intergenic_exp <- read.table("inputs/backmapping_stats/intergenic_allTissues_isoformTPM", header=T, stringsAsFactors=F)
tissue_specific_exp <- read.table("inputs/backmapping_stats/allTissues_isoformTPM", header=T, stringsAsFactors=F)
#remove mt entries
tissue_specific_exp <- tissue_specific_exp[-c(1,2),]
#combine
tissue_specific_all_exp <- rbind(tissue_specific_intergenic_exp,tissue_specific_exp)
refined_nolncRNA_exp <- merge(refined_nolncRNA,tissue_specific_all_exp, by.x="V4",by.y="isoformName")
refined_nolncRNA_exp <- refined_nolncRNA_exp[ , c("V4","BrainStem", "Cerebellum",  "Embryo.ICM", "Embryo.TE",  "Muscle",  "Retina",  "Skin",  "SpinalCord")] 

#attach tissue specific expression to our lncRNA
lncRNA_exp <- merge(all,tissue_specific_all_exp, by.x="V4",by.y="isoformName")
lncRNA_exp <- lncRNA_exp[ , c("V4","BrainStem", "Cerebellum",  "Embryo.ICM", "Embryo.TE",  "Muscle",  "Retina",	"Skin",	"SpinalCord")] 
#combine all
total_exp <- rbind(data.frame(id="lncRNA",lncRNA_exp),
                   data.frame(id="annotated",refined_nolncRNA_exp))
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

#subset lncRNA from each tissue
BS_lncRNA<-subset(BrainStem,id %in% "lncRNA")
C_lncRNA<-subset(Cerebellum,id %in% "lncRNA")
EICM_lncRNA<-subset(Embryo.ICM,id %in% "lncRNA")
ETE_lncRNA<-subset(Embryo.TE,id %in% "lncRNA")
M_lncRNA<-subset(Muscle,id %in% "lncRNA")
R_lncRNA<-subset(Retina,id %in% "lncRNA")
S_lncRNA<-subset(Skin,id %in% "lncRNA")
SC_lncRNA<-subset(SpinalCord,id %in% "lncRNA")

#make a cutoff of 0.1 TPM exp to be considered expressed in tissue
BS_lncRNA_cut<-subset(BS_lncRNA,BrainStem > 0.1)
#rownames(BS_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_BS<-sum(BS_lncRNA_cut[["BrainStem"]])
C_lncRNA_cut<-subset(C_lncRNA,Cerebellum > 0.1)
#rownames(C_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_C<-sum(C_lncRNA_cut[["Cerebellum"]])
EICM_lncRNA_cut<-subset(EICM_lncRNA,Embryo.ICM > 0.1)#
#rownames(EICM_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_EICM<-sum(EICM_lncRNA_cut[["Embryo.ICM"]])
ETE_lncRNA_cut<-subset(ETE_lncRNA,Embryo.TE > 0.1)
#rownames(ETE_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_ETE<-sum(ETE_lncRNA_cut[["Embryo.TE"]])
M_lncRNA_cut<-subset(M_lncRNA,Muscle > 0.1)
#rownames(M_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_M<-sum(M_lncRNA_cut[["Muscle"]])
R_lncRNA_cut<-subset(R_lncRNA,Retina > 0.1)
#rownames(R_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_R<-sum(R_lncRNA_cut[["Retina"]])
S_lncRNA_cut<-subset(S_lncRNA,Skin > 0.1)
#rownames(S_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_S<-sum(S_lncRNA_cut[["Skin"]])
SC_lncRNA_cut<-subset(SC_lncRNA,SpinalCord > 0.1)
#rownames(SC_lncRNA_cut) <- c()
cumulative_lncRNA_TPM_SC<-sum(SC_lncRNA_cut[["SpinalCord"]])

#now do the same for the genes
#subset PCGs from each tissue
BS_genes<-subset(BrainStem,id %in% "annotated")
C_genes<-subset(Cerebellum,id %in% "annotated")
EICM_genes<-subset(Embryo.ICM,id %in% "annotated")
ETE_genes<-subset(Embryo.TE,id %in% "annotated")
M_genes<-subset(Muscle,id %in% "annotated")
R_genes<-subset(Retina,id %in% "annotated")
S_genes<-subset(Skin,id %in% "annotated")
SC_genes<-subset(SpinalCord,id %in% "annotated")

#make a cutoff of 0.1 TPM exp to be considered expressed in tissue
BS_genes_cut<-subset(BS_genes,BrainStem > 0.1)
#rownames(BS_genes_cut) <- c()
cumulative_genes_TPM_BS<-sum(BS_genes_cut[["BrainStem"]])
C_genes_cut<-subset(C_genes,Cerebellum > 0.1)
#rownames(C_genes_cut) <- c()
cumulative_genes_TPM_C<-sum(C_genes_cut[["Cerebellum"]])
EICM_genes_cut<-subset(EICM_genes,Embryo.ICM > 0.1)
#rownames(EICM_genes_cut) <- c()
cumulative_genes_TPM_EICM<-sum(EICM_genes_cut[["Embryo.ICM"]])
ETE_genes_cut<-subset(ETE_genes,Embryo.TE > 0.1)
#rownames(ETE_genes_cut) <- c()
cumulative_genes_TPM_ETE<-sum(ETE_genes_cut[["Embryo.TE"]])
M_genes_cut<-subset(M_genes,Muscle > 0.1)
#rownames(M_genes_cut) <- c()
cumulative_genes_TPM_M<-sum(M_genes_cut[["Muscle"]])
R_genes_cut<-subset(R_genes,Retina > 0.1)
#rownames(R_genes_cut) <- c()
cumulative_genes_TPM_R<-sum(R_genes_cut[["Retina"]])
S_genes_cut<-subset(S_genes,Skin > 0.1)
#rownames(S_genes_cut) <- c()
cumulative_genes_TPM_S<-sum(S_genes_cut[["Skin"]])
SC_genes_cut<-subset(SC_genes,SpinalCord > 0.1)
#rownames(SC_genes_cut) <- c()
cumulative_genes_TPM_SC<-sum(SC_genes_cut[["SpinalCord"]])

#making lists for pie
BrainStem = list(c(cumulative_lncRNA_TPM_BS,cumulative_genes_TPM_BS))
Cerebellum = list(c(cumulative_lncRNA_TPM_C,cumulative_genes_TPM_C))
Embryo.ICM = list(c(cumulative_lncRNA_TPM_EICM,cumulative_genes_TPM_EICM))
Embryo.TE = list(c(cumulative_lncRNA_TPM_ETE,cumulative_genes_TPM_ETE))
Muscle = list(c(cumulative_lncRNA_TPM_M,cumulative_genes_TPM_M))
Retina = list(c(cumulative_lncRNA_TPM_R,cumulative_genes_TPM_R))
Skin = list(c(cumulative_lncRNA_TPM_S,cumulative_genes_TPM_S))
SpinalCord = list(c(cumulative_lncRNA_TPM_SC,cumulative_genes_TPM_SC))


#Pies only including genes and lncRNA
library(caroline)
pdf("Fig4A.pdf")
par(lwd = 2.5)
pies(
  list(
    Brainstem=nv(c(634533,752408),c('lncRNA','annotated')),
    Cerebellum=nv(c(962495,743612),c('lncRNA','annotated')),
    Embryo_ICM=nv(c(637339,577346),c('lncRNA','annotated')),
    Embryo_TE=nv(c(586020,609540),c('lncRNA','annotated')),
    Muscle=nv(c(118968,848670),c('lncRNA','annotated')),
    Retina=nv(c(771653,933294),c('lncRNA','annotated')),
    Skin=nv(c(466649,941957),c('lncRNA','annotated')),
    Spinal_cord=nv(c(779894,792801),c('lncRNA','annotated'))),
  x0=c(57892,58829,51222,47204,40102,46528,45107,58897),
  y0=c(18660,18917,16715,20800,10540,12942,13749,18567),
  radii=4, border=c('purple','purple','black','black','blue','blue','blue','purple'))#,show.labels=T) # to see labels add ",show.labels=T"
dev.off()

#B) tissue-specific heatmap
#Melt data for manipulation
library(reshape2)
melted <- melt(lncRNA_exp, id.vars=c("V4"))
#Calculate the sum(TPM) and STDEV of each gene per tissue
melted_new<- ddply(melted, c("V4"), summarise,
                   sum = sum(value), sd = sd(value),
                   sem = sd(value)/sqrt(length(value)))
#Add the column of sum and sd to the original TPM values table
complete <- merge(lncRNA_exp,melted_new,by="V4")
#Subset data based on if sum>50 and sd>50
sub <- subset(complete, c(sum > 100 & sd > 50))
#make row.names the geneName
rownames(sub)<-sub$V4
#rownames(lncRNA_exp)<-lncRNA_exp$V4
# making the matrix for the heatmap
#disable scientific notation so no "e+/-"
options("scipen"=100, "digits"=4)
#datanumbers <- data.matrix(lncRNA_exp[,2:9])
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
par(mar=c(7,4,4,2)+0.1) 
png(filename='heatmap_manual_lncRNA.png', width=800, height=750)
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
hmap_order <- data.frame(datanumbers_smalls[rev(hr$labels[hr$order]), hc$labels[hc$order]])
rownames(hmap_order) -> hmap_order$V4
#merge with lncRNA bed
lncRNA_bed <- read.table("lncRNA_final.bed", header=F,stringsAsFactors=FALSE)

##using a merge function provided by Tal Galili (thank you!) to retain order
merge.with.order <- function(x,y, ..., sort = T, keep_order)
{
  add.id.column.to.data <- function(DATA)
  {
    data.frame(DATA, id... = seq_len(nrow(DATA)))
  }
  order.by.id...and.remove.it <- function(DATA)
  {
    if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")  
    ss_r <- order(DATA$id...)
    ss_c <- colnames(DATA) != "id..."
    DATA[ss_r, ss_c]
  }
  if(!missing(keep_order))
  {
    if(keep_order == 1) return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
    if(keep_order == 2) return(order.by.id...and.remove.it(merge(x=x,y=add.id.column.to.data(y),..., sort = FALSE)))
    warning("The function merge.with.order only accepts NULL/1/2 values for the keep_order variable")
  } else {return(merge(x=x,y=y,..., sort = sort))}
}
#
hmap_order_coord <- merge.with.order( hmap_order, lncRNA_bed, by='V4', sort=F ,keep_order = 1)
write.csv(datanumbers_smalls[rev(hr$labels[hr$order]), hc$labels[hc$order]],"tissue_4_hmap_order.csv")
write.csv(hmap_order_coord ,"tissue_lncRNA_ordered_coord.csv")

hmap_order_coord <- hmap_order_coord[ ,c(1:10)]
#elt data so I can use group by
melted_hmap_order_coord <- melt(hmap_order_coord, id.vars=c("V4","V1"))
#to find which chromosome is most frequent
freq_chr2 <- as.table(with(melted_hmap_order_coord,by(V1,variable,function(xx)names(which.max(table(V1))))))

#C) unique vs absent
library(ggplot2)
###plotting figures with varying threshold for absent vs unique lncRNA
data_0.1<-read.table("tissueSpecificSummary_cutoff.0.1")
data_changed_0.1 <- cbind(as.data.frame(data_0.1[1:4,]),as.data.frame(data_0.1[5:8,]),as.data.frame(data_0.1[9:12,]),as.data.frame(data_0.1[13:16,]),as.data.frame(data_0.1[17:20,]),as.data.frame(data_0.1[21:24,]),as.data.frame(data_0.1[25:28,]),as.data.frame(data_0.1[29:32,]))   
data_changed_0.1 <- sapply(data_changed_0.1, as.character)
colnames(data_changed_0.1) <- data_changed_0.1[1,]
data_changed_0.1 <- as.data.frame(data_changed_0.1[-1,])
data_0.1_table <-as.data.frame(t(data_changed_0.1))
rownames(data_changed_0.1) <- c("total","unique_lncRNA","not_unique_lncRNA")
data_changed_0.1 <- data_changed_0.1[c(2,3), ]
data <-as.data.frame(t(data_changed_0.1))
write.table(data, "tissue_lncRNA_1.txt")
data_01 <- read.table("tissue_lncRNA_1.txt",stringsAsFactors=FALSE)
data_01$not_unique_lncRNA <- data_01$not_unique_lncRNA * -1
Absent_U_data_0.1 <- as.data.frame(data_01$not_unique_lncRNA)
rownames(Absent_U_data_0.1) <- rownames(data_01) 
U_data_1 <- as.data.frame(data_01$unique_lncRNA)
rownames(U_data_1) <- rownames(data_01)

# to get cumulative TPM of those expressed
BrainStem <- read.table("uniqExp_lncRNA/BrainStem.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE)
BrainStem <- BrainStem[2]
BrainStem <- sum(BrainStem)
Cerebellum <- read.table("uniqExp_lncRNA/Cerebellum.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
Cerebellum <- Cerebellum[3]
Cerebellum <- sum(Cerebellum)
Embryo.ICM <- read.table("uniqExp_lncRNA/Embryo.ICM.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
Embryo.ICM <- Embryo.ICM[4]
Embryo.ICM <- sum(Embryo.ICM)
Embryo.TE <- read.table("uniqExp_lncRNA/Embryo.TE.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
Embryo.TE <- Embryo.TE[5]
Embryo.TE <- sum(Embryo.TE)
Muscle <- read.table("uniqExp_lncRNA/Muscle.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
Muscle <- Muscle[6]
Muscle <- sum(Muscle)
Retina <- read.table("uniqExp_lncRNA/Retina.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
Retina <- Retina[7]
Retina <- sum(Retina)
Skin <- read.table("uniqExp_lncRNA/Skin.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
Skin <- Skin[8]
Skin <- sum(Skin)
SpinalCord <- read.table("uniqExp_lncRNA/SpinalCord.isoform.expressed_uniqely_cutoff.0.1",stringsAsFactors=FALSE) 
SpinalCord <- SpinalCord[9]
SpinalCord <- sum(SpinalCord)
sums <- data.frame(c(BrainStem,Cerebellum,Embryo.ICM,Embryo.TE,Muscle,Retina,Skin,SpinalCord))
row.names(sums) <- c("BrainStem","Cerebellum","Embryo.ICM","Embryo.TE","Muscle","Retina","Skin","SpinalCord")
names(sums) <- c("sum")

png(filename='Fig4C.png', width=800, height=750)
ggplot() +
  geom_bar(data=U_data_1, aes(x=rownames(data_01),y=data_01$unique_lncRNA,color="aliceblue"), stat="identity") +
  geom_bar(data=Absent_U_data_0.1, aes(x=rownames(data_01),y=data_01$not_unique_lncRNA,color="red"),stat = "identity") + 
  ylab("Number of lncRNA") + scale_color_discrete(name="Unique lncRNA",
                                                  labels=c("present","absent")) + xlab("Tissue") +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 12)) +
  theme(axis.text.x = element_text(colour="black", size = 9)) +
  theme(axis.title = element_text(colour="black", size = 14)) +
  geom_line(data=sums, aes(x=row.names(sums),y=sum / 5, group=1),colour="green")
graphics.off()
