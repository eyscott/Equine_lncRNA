##making a distribution plot with transcripts and TPM
setwd("~/Desktop/lncRNA")
total_exp <- read.table("total_tissue_all_TPM", header=T, stringsAsFactors=F)
setwd("~/Dropbox/Horse_Transcriptome/downloads")
refined_T_exp <- read.table("backmapping_stats/allTissues_isoformTPM", header=T, stringsAsFactors=F)
refined_G_exp <- read.table("backmapping_stats/allTissues_geneTPM", header=T, stringsAsFactors=F)
refined_T_exp <- refined_T_exp[-c(1,2),]
refined_G_exp <- refined_G_exp[-c(1,2),]
rownames(refined_T_exp) <- c()
rownames(refined_G_exp) <- c()
#Making for known lncRNA
setwd("~/Desktop/lncRNA")
known_and_unfiltered <- read.table("known_lncRNA_exp.bed", header=T, stringsAsFactors=F)#15584
names(known_and_unfiltered)[1] <- c("transcriptName")
#also for our identified lncRNA
setwd("~/Desktop/lncRNA")
our_lncRNA <- read.table("lncRNA_exp.txt", header=T, stringsAsFactors=F)#18025
#combine ours and known lncRNA
all_lncRNA <- rbind(known_and_unfiltered,our_lncRNA)
rownames(all_lncRNA) <- c()
all_lncRNA <- all_lncRNA[!duplicated(all_lncRNA),]
#without muscle and retina and skin
all_lncRNA_noMorSorR <-data.frame(all_lncRNA[ ,-c(6:9)])
#reshaping data for more easy manipulation
library(reshape2)
#total_exp <- total_exp[-1]
total_melt <- melt(total_exp, id.vars=c("transcriptName"))
refined_melt <- melt(refined_T_exp, id.vars=c("isoformName"))
refined_G_melt <- melt(refined_G_exp, id.vars=c("geneName"))
known_and_unfiltered_melt <- melt(known_and_unfiltered, id.vars=c("transcriptName"))
our_lncRNA_melt <- melt(our_lncRNA, id.vars=c("transcriptName"))
all_lncRNA_melt <- melt(all_lncRNA, id.vars=c("transcriptName"))
all_lncRNA_noMorSorR_melt <- melt(all_lncRNA_noMorSorR, id.vars=c("transcriptName"))
#variance stabilizing data
library(plyr)
melted_new<- ddply(total_melt, c("transcriptName","variable"), summarise,
                   stable=log10(value + 1))
refined_melted_new<- ddply(refined_melt, c("isoformName","variable"), summarise,
                   stable=log10(value + 1))
refined_G_melted_new<- ddply(refined_G_melt, c("geneName","variable"), summarise,
                           stable=log10(value + 1))
known_and_unfiltered_melt_new<- ddply(known_and_unfiltered_melt, c("transcriptName","variable"), summarise,
                                      stable=log10(value + 1))
our_lncRNA_melt_new<- ddply(our_lncRNA_melt, c("transcriptName","variable"), summarise,
                                      stable=log10(value + 1))
all_lncRNA_melt_new<- ddply(all_lncRNA_melt, c("transcriptName","variable"), summarise,
                            stable=log10(value + 1))
all_lncRNA_noMorSorR_melt_new<- ddply(all_lncRNA_noMorSorR_melt, c("transcriptName","variable"), summarise,
                            stable=log10(value + 1))

#Make the density plot for each
library(RColorBrewer)
library(ggplot2)
my.cols <- brewer.pal(8, "Set1")
ggplot(melted_new, aes(x=stable)) + geom_density(aes(group=variable,colour=variable)) +
  xlim(-1,3) + xlab("log2(TPM+1)") + scale_colour_manual(values = my.cols)
#for just refined annotated transcripts
r <- ggplot(refined_melted_new, aes(x=stable)) + geom_density(aes(group=variable,colour=variable)) +
   xlab("log10(TPM+1)") + scale_colour_manual(values = my.cols) +    
  xlim(-1,3)

t <- ggplot(refined_melt, aes(x=value)) + geom_density(aes(group=variable,colour=variable)) +
  xlab("TPM") + scale_colour_manual(values = my.cols)  
#xlim(0,10)

#for just refined annotated genes
g <- ggplot(refined_G_melted_new, aes(x=stable)) + geom_density(aes(group=variable,colour=variable)) +
  xlab("log2(TPM+1)") + scale_colour_manual(values = my.cols) + xlim(-1,10)

#for known lncRNA
ggplot(known_and_unfiltered_melt, aes(x=value)) + geom_density(aes(group=variable,colour=variable)) +
  xlab("TPM") + xlim(-1,3) + scale_colour_manual(values = my.cols) 
#variance stabilized version
l <- ggplot(known_and_unfiltered_melt_new, aes(x=stable)) + geom_density(aes(group=variable,colour=variable)) +
  xlab("log2(TPM+1)") + scale_colour_manual(values = my.cols) 
  # xlim(0,10)

#for our lncRNA
ggplot(our_lncRNA_melt, aes(x=value)) + geom_density(aes(group=variable,colour=variable)) +
  xlab("TPM") + xlim(-1,3) + scale_colour_manual(values = my.cols) 
#variance stabilized version
o <- ggplot(our_lncRNA_melt_new, aes(x=stable)) + geom_density(aes(group=variable,colour=variable)) +
  xlab("log10(TPM+1)") + xlim(-1,3) + scale_colour_manual(values = my.cols) 
#for our lncRNA and known
ggplot(all_lncRNA_melt, aes(x=value)) + geom_density(aes(group=variable,colour=variable)) +
  xlab("TPM") + xlim(-1,3) + scale_colour_manual(values = my.cols) 
#variance stabilized version
a <- ggplot(all_lncRNA_melt_new, aes(x=stable)) + geom_density(aes(group=variable,colour=variable)) +
  xlab("log10(TPM+1)") + scale_colour_manual(values = my.cols) + 
   xlim(-0.5,1)

#muscle and retina kind of mess this up bc of scale
#without muscle and retina
x <- ggplot(all_lncRNA_noMorSorR_melt_new, aes(x=stable)) + geom_density(aes(group=variable,colour=variable)) +
  xlab("log10(TPM+1)") + xlim(-1,3) + scale_colour_manual(values = my.cols) 

#getting all the density values
refined_points<-print(r)
refined_points<-refined_points$data[[1]]
write.table(refined_points,"refined_points.txt")

lncRNA_points<-print(l)
lncRNA_points<-lncRNA_points$data[[1]]

all_lncRNA_points<-print(a)
all_lncRNA_points<-all_lncRNA_points$data[[1]]



#make a function for calculating mode with colour as a factor
#apparently while subsetting you can not have a mix of positive anf nevatige numbers
colours <- factor(refined_points$colour)
##377EB8 #4DAF4A #984EA3 #A65628 #E41A1C #F781BF #FF7F00 ##FFFF33
#rename colour factors to tissue names
levels(colours) <- list(BrainStem="#377EB8", Cerebellum="#4DAF4A", Embryo.ICM="#984EA3", Embryo.TE="#A65628", Muscle="#E41A1C", Retina="#F781BF", Skin="#FF7F00", SpinalCord="#FFFF33")

#obtaining the TPM that has highest density for each colour (=tissue)
getmode <- function(v) {
  v[which.max(v$density),]
}
result_refined <- by(refined_points,colours, getmode, simplify=T)
print(result_refined)

result_lncRNA <- by(lncRNA_points,colours, getmode, simplify=T)
print(result_lncRNA)

result_all_lncRNA <- by(all_lncRNA_points,colours, getmode, simplify=T)
print(result_all_lncRNA)

#Obtain interval datapoints to make P(detection) curve
#calculate areas with 0.1 TPM intervals
getAreas_please <- function(v) {for (i in 1:dim(v)) {
  I<-v$x[i]-v$x[i-1]
  A<-v$density*I
}
v$area<-print(A)
}
#these are the areas for each TPM
intervals_refined <- getAreas_please(refined_points)
intervals_lncRNA <- getAreas_please(lncRNA_points)
intervals_all_lncRNA <- getAreas_please(all_lncRNA_points)
#now attach these back to the TPM values and tissue categories
refined_area<-data.frame(cbind(refined_points$x,intervals_refined,refined_points$colour))
lncRNA_area<-data.frame(cbind(lncRNA_points$x,intervals_lncRNA,lncRNA_points$colour))
all_lncRNA_area<-data.frame(cbind(all_lncRNA_points$x,intervals_all_lncRNA,all_lncRNA_points$colour))
#make sure all values besides factors are numeric so we can perform calculations
refined_area[, 1] <- as.numeric(as.character( refined_area[, 1] ))
refined_area[, 2] <- as.numeric(as.character( refined_area[, 2] ))
lncRNA_area[, 1] <- as.numeric(as.character( lncRNA_area[, 1] ))
lncRNA_area[, 2] <- as.numeric(as.character( lncRNA_area[, 2] ))
all_lncRNA_area[, 1] <- as.numeric(as.character( all_lncRNA_area[, 1] ))
all_lncRNA_area[, 2] <- as.numeric(as.character( all_lncRNA_area[, 2] ))

#divide these intervals by AUC
#sum all the areas for each tissue, so again set colour(=tissue) as factor
colours <- factor(refined_area$V3)
#rename colour factors to tissue names
levels(colours) <- list(BrainStem="#377EB8", Cerebellum="#4DAF4A", Embryo.ICM="#984EA3", Embryo.TE="#A65628", Muscle="#E41A1C", Retina="#F781BF", Skin="#FF7F00", SpinalCord="#FFFF33")
#calculate the AUC
AUC_refined <- aggregate(intervals_refined ~ colours, refined_area, sum)
AUC_lncRNA <- aggregate(intervals_lncRNA ~ colours, lncRNA_area, sum)
AUC_all_lncRNA <- aggregate(intervals_all_lncRNA ~ colours, all_lncRNA_area, sum)
#divide each interval by the respective AUC
refined <- split(refined_area, refined_area$V3)
str(refined)
# make each into a data.frame
Y <- lapply(seq_along(refined), function(x) as.data.frame(refined[[x]])) 
#Assign the dataframes in the list Y to individual objects and divide them by their
#respective AUC
library(dplyr)
#get the AUC values from above and apply manually
BS_r <- Y[[1]]
BS_r_f <- mutate(BS_r, P = intervals_refined / 0.9483482)
C_r <- Y[[2]]
C_r_f <- mutate(C_r, P = intervals_refined / 0.9009495)
EICM_r <- Y[[3]]
EICM_r_f <- mutate(EICM_r, P = intervals_refined / 0.8781400)
ETE_r <- Y[[4]]
ETE_r_f <- mutate(ETE_r, P = intervals_refined / 0.8627394)
M_r <- Y[[5]]
M_r_f <- mutate(M_r, P = intervals_refined / 0.9398763)
R_r <- Y[[6]]
R_r_f <- mutate(R_r, P = intervals_refined / 0.9473616)
S_r <- Y[[7]]
S_r_f <- mutate(S_r, P = intervals_refined / 0.8574869)
SC_r <- Y[[8]]
SC_r_f <- mutate(SC_r, P = intervals_refined / 0.8524720)
P_refined <- rbind(BS_r_f,C_r_f,EICM_r_f,ETE_r_f,M_r_f,R_r_f,S_r_f,SC_r_f)

#do for known lncRNA
lncRNA <- split(lncRNA_area, lncRNA_area$V3)
str(lncRNA)
# make each into a data.frame
X <- lapply(seq_along(lncRNA), function(x) as.data.frame(lncRNA[[x]])) 
#Assign the dataframes in the list Y to individual objects
BS_l <- X[[1]]
BS_l_f <- mutate(BS_l, P = intervals_lncRNA / 0.9553902)
C_l <- X[[2]]
C_l_f <- mutate(C_l, P = intervals_lncRNA / 0.8614871)
EICM_l <- X[[3]]
EICM_l_f <- mutate(EICM_l, P = intervals_lncRNA / 0.8403547)
ETE_l <- X[[4]]
ETE_l_f <- mutate(ETE_l, P = intervals_lncRNA / 0.8466310)
M_l <- X[[5]]
M_l_f <- mutate(M_l, P = intervals_lncRNA / 0.9484222)
R_l <- X[[6]]
R_l_f <- mutate(R_l, P = intervals_lncRNA / 0.9602346)
S_l <- X[[7]]
S_l_f <- mutate(S_l, P = intervals_lncRNA / 0.8331962)
SC_l <- X[[8]]
SC_l_f <- mutate(SC_l, P = intervals_lncRNA / 0.8763995)
P_lncRNA <- rbind(BS_l_f,C_l_f,EICM_l_f,ETE_l_f,M_l_f,R_l_f,S_l_f,SC_l_f)

#now do for all lncRNA
#do for known lncRNA
all_lncRNA <- split(all_lncRNA_area, all_lncRNA_area$V3)
str(all_lncRNA)
# make each into a data.frame
Z <- lapply(seq_along(all_lncRNA), function(x) as.data.frame(all_lncRNA[[x]])) 
#Assign the dataframes in the list Y to individual objects
BS_la <- Z[[1]]
BS_la_f <- mutate(BS_la, P = intervals_all_lncRNA / 0.9397867)
C_la <- Z[[2]]
C_la_f <- mutate(C_la, P = intervals_all_lncRNA / 0.9477055)
EICM_la <- Z[[3]]
EICM_la_f <- mutate(EICM_la, P = intervals_all_lncRNA / 0.8520282)
ETE_la <- Z[[4]]
ETE_la_f <- mutate(ETE_la, P = intervals_all_lncRNA / 0.7612786)
M_la <- Z[[5]]
M_la_f <- mutate(M_la, P = intervals_all_lncRNA / 0.8951038)
R_la <- Z[[6]]
R_la_f <- mutate(R_la, P = intervals_all_lncRNA / 0.8841801)
S_la <- Z[[7]]
S_la_f <- mutate(S_la, P = intervals_all_lncRNA / 1.6235155)
SC_la <- Z[[8]]
SC_la_f <- mutate(SC_la, P = intervals_all_lncRNA / 0.7387983)
P_all_lncRNA <- rbind(BS_la_f,C_la_f,EICM_la_f,ETE_la_f,M_la_f,R_la_f,S_la_f,SC_la_f)

#make a scatter plot with fit line of these fractions against TPM 
#v1=TPM=x-axis, P=probability=y
colours <- factor(P_refined$V3)
##377EB8 #4DAF4A #984EA3 #A65628 #E41A1C #F781BF #FF7F00 ##FFFF33
#rename colour factors to tissue names
levels(colours) <- list(BrainStem="#377EB8", Cerebellum="#4DAF4A", Embryo.ICM="#984EA3", Embryo.TE="#A65628", Muscle="#E41A1C", Retina="#F781BF", Skin="#FF7F00", SpinalCord="#FFFF33")

library(ggplot2)
#define xlim() by the mode/ highest peak in density plot
r2<- ggplot(P_refined,aes(V1,P,colour=V3)) +
  geom_line() + xlab("log10(TPM+1)") + ylab("Probability of Detection") +
  scale_color_discrete(name="Tissues",
                       labels=c("BrainStem", "Cerebellum",  "Embryo.ICM", "Embryo.TE",	"Muscle",	"Retina",	"Skin",	"SpinalCord")) 
 

l<- ggplot(P_lncRNA,aes(V1,P,colour=V3)) +
  geom_line() + xlim(0,.035) + xlab("log2(TPM+1)") + ylab("Probability of Detection") +
  scale_color_discrete(name="Tissues",
                       labels=c("BrainStem", "Cerebellum",  "Embryo.ICM", "Embryo.TE",  "Muscle",	"Retina",	"Skin",	"SpinalCord"))


la<- ggplot(P_all_lncRNA,aes(V1,P,colour=V3)) +
  geom_line() + xlim(0,2) + ylim(0,0.3) +
  xlab("log2(TPM+1)") + ylab("Probability of Detection") +
  scale_color_discrete(name="Tissues",
                       labels=c("BrainStem", "Cerebellum",  "Embryo.ICM", "Embryo.TE",  "Muscle",	"Retina",	"Skin",	"SpinalCord"))

#find the max P(detection) for the mode of refined genes
#and lncRNA for each tissue
getTPM <- function(v) {
  v[which.max(v$P),]
}

library(dplyr)
P_refined_density <- cbind(P_refined, refined_freq)
refined_P_max <- data.frame(do.call("rbind", by(P_refined,colours, getTPM, simplify=T)))

lncRNA_P_max <- data.frame(do.call("rbind", by(P_lncRNA,colours, getTPM, simplify=T)))

all_lncRNA_P_max <- data.frame(do.call("rbind", by(P_all_lncRNA,colours, getTPM, simplify=T)))

P_detect_refined_lncRNA <- cbind(refined_P_max,all_lncRNA_P_max)
P_detect_refined_lncRNA <- P_detect_refined_lncRNA[ ,c(1,4,5,8)]
names(P_detect_refined_lncRNA) <- c("TPM_refined_log","P_refined","TPM_lncRNA_log","P_lncRNA")

setwd("~/Desktop/lncRNA")
write.table(P_detect_refined_lncRNA,"P_detection_table.txt")

#plotting relationship between P(detecting genes) vs P(detecting lncRNA)
setwd("~/Desktop/lncRNA")
P_detections <- read.table("P_detection_table.txt", header=T, stringsAsFactors=F)
rownames(P_detections)->P_detections$Tissue
library(RColorBrewer)
library(ggplot2)
my.cols <- brewer.pal(8, "Set1")
d <- ggplot(P_detections) +
  geom_point( aes(x=P_refined,y=P_lncRNA,colour=tissue)) + 
  scale_colour_manual(values = my.cols) + 
  xlab("P(detecting mode expression of PCG)") +
  ylab("log10(P(detecting mode expression of lncRNA))") +
  scale_y_log10(breaks = c(0.01,0.02,0.03,0.04,0.05,0.06,0.1,0.2,0.3,0.4,0.5,0.6))
