setwd("C:/Users/aki/Desktop/S.latissima/analysis/")

require(ggplot2)

pop<-c(rep("AKR",7),rep("GAR",8),rep("GRO",8),rep("GRU",10),rep("HEL",8),rep("POR",8))

z<-read.csv("seqsummary.csv")

ggplot(z, aes(x="", y=Avg_seq_site, fill=Site)) +
  geom_bar(stat="identity", width=1) +
  ggtitle("Average number of sequence reads")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round((Avg_seq_site/sum(Avg_seq_site)*100),1), "%")),position = position_stack(vjust=0.5),size=3.5) +
  labs(x = NULL, y = NULL, fill = NULL)+
  theme_void()

ggplot(z, aes(x="", y=Avg_base_site, fill=Site)) +
  geom_bar(stat="identity", width=1) +
  ggtitle("Average number of sequenced bases")+
  coord_polar("y", start=0)+
  geom_text(aes(label = paste0(round((Avg_base_site/sum(Avg_base_site)*100),1), "%"),x = 1.65),position = position_stack(vjust=0.5),size=3.5) +
  labs(x = NULL, y = NULL, fill = NULL)+
  theme_void()

x<-read.csv("alignments.csv")
x$pop<-c(rep("AKR",7),rep("GAR",8),rep("GRO",8),rep("GRU",10),rep("HEL",8),rep("POR",8))
x<-x[order(x$Reads_aligned,decreasing = T),]
x$Sample<-factor(x$Sample,levels = x$Sample)

ggplot(data=x)+
  geom_col(aes(x = Sample, y = Reads_aligned))+
  ylab("Number of aligned reads")+
  ggtitle("Reads aligned prior to duplicate removal")+
  ylim(0,220000000)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.25))

y<-x[order(x$Reads_aligned_postMD,decreasing = T),]
y$Sample<-factor(y$Sample,levels = y$Sample)

ggplot(data=y)+
  geom_col(aes(x = Sample, y = Reads_aligned_postMD))+
  ylab("Number of aligned reads")+
  ggtitle("Reads aligned post duplicate removal")+
  ylim(0,220000000)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.25))

round((aggregate(x$Reads_aligned_postMD,list(x$pop),FUN=sum)[,2]/sum(x$Reads_aligned_postMD)*100),2)
ggplot(x, aes(x="", y=Reads_aligned_postMD, fill=pop,group=pop)) +
  geom_bar(stat="identity", width=1) +
  #ggtitle("Average number of sequence reads")+
  coord_polar("y", start=0)+
  #stat_summary(aes(label=..y..), fun.y=mean, geom="text", size=8)+
  #geom_text(aes(label = paste0(),position = position_stack(vjust=0.5),size=3.5) +
  labs(x = NULL, y = NULL, fill = NULL)+
  theme_void()

#rm(list=ls())

indmiss<-read.table("out.imiss",sep="\t",header = T)
indmiss$INDV

ind<-c(paste("AKR",c(1,2,4:8), sep=""),paste("GAR",c(2:9), sep=""),paste("GRO",c(1:7,9), sep=""),paste("GRU",c(1,10,2:9), sep=""),paste("HEL",c(2:9), sep=""),paste("POR",c(1:4,6:9), sep=""))

indmiss$Sample<-ind
indmiss<-indmiss[order(indmiss$F_MISS,decreasing = T),]
indmiss$Sample<-factor(indmiss$Sample,levels = indmiss$Sample)
ggplot(data=indmiss[indmiss$F_MISS<1,])+
  geom_col(aes(x = Sample, y = (F_MISS)))+
  ylab("Fraction of missing SNPs")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = -45, vjust = 0.25))

require(vcfR)
#x<-read.vcfR("S.laminaria_bowtie_samtools.mQ30mMQ40.vcf.gz")
x<-read.vcfR("S.latissima_bowtie_samtools.mQ30mMQ40.g5mac3dp3maf05DP20p5b.recode.vcf")
geno<-extract.gt(x)
position <- getPOS(x) # Positions in bp
chromosome <- getCHROM(x) # Chromosome information
#colnames(geno)<-c(paste("AKR",c(1,2,4:8), sep=""),paste("GAR",c(2:9), sep=""),paste("GRO",c(1:7,9), sep=""),paste("GRU",c(1,10,2:9), sep=""),paste("HEL",c(2:9), sep=""),paste("POR",c(1:4,6:9), sep=""))
#gt<-data.frame(geno)

#write.table(noquote(cbind(colnames(geno),pop)),file="pop.txt",sep=" ",col.names = F,row.names = F,quote = F)

sum(is.na(geno))/(nrow(geno)*ncol(geno))
summary(as.factor(geno))

genoV<-NULL
snp_name<-NULL
for(i in 1:nrow(geno)){
  if(!"0/2"%in%geno[i,]){
  #if(sum(is.na(geno[i,]))==0){
    #if(length(unique(geno[i,]))!=1){
      genoV<-rbind(genoV,geno[i,])
      snp_name<-c(snp_name,row.names(geno)[i])
    #}
  }
}

badloci<-NULL
for(i in 1:nrow(geno)){
  if("0/2"%in%geno[i,]){
    badloci<-c(badloci,row.names(geno)[i])
  }
}

write.table(badloci,file="weirdloci",sep=" ",col.names = F,row.names = F,quote = F)

summary(as.factor(genoV))
colnames(genoV)<-colnames(geno)
row.names(genoV)<-snp_name

write.table(genoV,file = "geno_clean.tsv", sep = "\t",row.names = T, col.names = T)

#require(BiocManager)
#BiocManager::install("qvalue")
#require(devtools)
#devtools::install_github("whitlock/OutFLANK")
require(OutFLANK)

#BiocManager::install("SNPRelate")
require(SNPRelate)

#remotes::install_github("privefl/bigsnpr")
#require(bigsnpr)

#install.packages("robust")
#require(robust)

#install.packages("bigstatsr")
#require(bigstatsr)

pos_loc <- row.names(genoV)

G <- matrix(NA, nrow = nrow(genoV), ncol = ncol(genoV),dimnames = list(pos_loc,colnames(genoV)) )

G[genoV %in% c("0/0")] <- 0
G[genoV %in% c("0/1", "1/0")] <- 1
G[genoV %in% c("1/1")] <- 2
G[genoV %in% NA] <- 9

SNPmat<-(t(G))

#modified function MakeDiploidFSTMat
MakeDiploidFSTMat<-function(SNPmat,locusNames,popNames){
  locusname <- unlist(locusNames)
  popname <- unlist(popNames)
  snplevs <- levels(as.factor(unlist(SNPmat)))
  if(any(!(snplevs%in%c(0,1,2,9)))==TRUE) {
    print("Error: Your snp matrix has a character other than 0,1,2 or 9")
    break
  }
  if (dim(SNPmat)[1] != length(popname)) {
    print("Error: your population names do not match your SNP matrix")
    break
  }
  if (dim(SNPmat)[2] != length(locusname)) {
    print("Error:  your locus names do not match your SNP matrix")
    break
  }
  writeLines("Calculating FSTs, may take a few minutes...")
  nloci <- length(locusname)
  FSTmat <- matrix(NA, nrow = nloci, ncol = 8)
  for (i in 1:nloci) {
    FSTmat[i, ] = unlist(getFSTs_diploids(popname, SNPmat[,i]))
    if (i%%10000 == 0) {
      print(paste(i, "done of", nloci))
    }
  }
  outTemp = as.data.frame(FSTmat)
  outTemp = cbind(locusname, outTemp)
  colnames(outTemp) = c("LocusName", "He", "FST", "T1", "T2", 
                        "FSTNoCorr", "T1NoCorr", "T2NoCorr", "meanAlleleFreq")
  return(outTemp)
}

getFSTs_diploids = function(popNameList, SNPDataColumn){  
  #eliminating the missing data for this locus
  popnames=unlist(as.character(popNameList))
  popNameTemp=popnames[which(SNPDataColumn!=9)]
  snpDataTemp=SNPDataColumn[SNPDataColumn!=9]
  
  HetCounts <- tapply(snpDataTemp, list(popNameTemp,snpDataTemp), length)
  HetCounts[is.na(HetCounts)] = 0
  
  #Case: all individuals are genetically identical at this locus
  if(dim(HetCounts)[2]==1){
    return (list(He=NA,FST=NA, T1=NA, T2=NA,FSTNoCorr=NA, T1NoCorr=NA, T2NoCorr=NA,meanAlleleFreq = NA))
  }
  
  if(dim(HetCounts)[2]==2){
    if(paste(colnames(HetCounts),collapse="")=="01"){HetCounts=cbind(HetCounts,"2"=0)}
    if(paste(colnames(HetCounts),collapse="")=="12"){HetCounts=cbind("0"=0,HetCounts)} 
    if(paste(colnames(HetCounts),collapse="")=="02"){HetCounts=cbind(HetCounts[,1],"1"=0, HetCounts[,2])}
  }
  
  out = WC_FST_Diploids_2Alleles(HetCounts)	
  return(out)
}

my_fst <- MakeDiploidFSTMat(SNPmat, locusNames = pos_loc, popNames = pop)

# If a locus has a much lower sample size compared to the rest, it could have a broader error distribution (and therefore incorrectly inferred as an outlier).
plot(my_fst$FST, my_fst$FSTNoCorr)
abline(a=0,b=1,col="red")

ind.keep<-1:nrow(G)

out_trim <- OutFLANK(my_fst[ind.keep,], NumberOfSamples=length(unique(pop)), qthreshold = 0.001, Hmin = 0.01)
str(out_trim)
head(out_trim$results)
summary(out_trim$results$OutlierFlag)
summary(out_trim$results$pvaluesRightTail)

OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
## Zoom in on right tail
OutFLANKResultsPlotter(out_trim , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.01, binwidth = 0.01, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)


#vcf.fn<-"S.laminaria_bowtie_samtools.mQ30mMQ40.vcf.gz"
#vcf.fn<-"S.laminaria_bowtie_samtools.mQ30mMQ40.g5mac3dp3indg95maf05p5b.recode.vcf"
vcf.fn<-"S.latissima_bowtie_samtools.mQ30mMQ40.g5mac3dp3maf05DP20p5b.recode.vcf"
#snpgdsVCF2GDS(vcf.fn, "S.laminaria_raw.gds",method = "copy.num.of.ref")
#snpgdsVCF2GDS(vcf.fn, "S.laminaria_fullfilter.gds",method = "copy.num.of.ref")
snpgdsVCF2GDS(vcf.fn, "S.laminaria_modfilter.gds",method = "copy.num.of.ref")

#Slam<-snpgdsOpen("S.laminaria_raw.gds")
#Slam<-snpgdsOpen("S.laminaria_fullfilter.gds")
Slam<-snpgdsOpen("S.laminaria_modfilter.gds")
#snpgdsClose(Slam)

pcaC <- snpgdsPCA(Slam,autosome.only=F)

pc.percent <- pcaC$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pcaC$sample.id,
                  EV1 = pcaC$eigenvect[,1],
                  EV2 = pcaC$eigenvect[,2],
                  EV3 = pcaC$eigenvect[,3],
                  EV4 = pcaC$eigenvect[,4],
                  #EV5 = pcaC$eigenvect[,5],
                  #EV6 = pcaC$eigenvect[,6],
                  stringsAsFactors = FALSE)
head(tab)

Per_exp<-head(round(pc.percent, 2))

pop<-c(rep("AKR",7),rep("GAR",8),rep("GRO",8),rep("GRU",10),rep("HEL",8),rep("POR",8))
#pop2<-c(rep("AKR",2),rep("GAR",4),rep("GRO",1),rep("HEL",4),rep("POR",7))

tab$pop<-pop
#tab$pop<-pop2
tab$ind<-ind #ind object from line 64
cl<-c("blue","cyan", "orange", "darkgreen","red", "black")
#cl<-c("blue","cyan", "orange", "darkgreen","red")
#okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#E15759")
shp<-c(0,1,2,3,4,5)
#shp<-c(0,1,2,3,4)
ggplot(tab, aes(EV1,EV2,color=pop,shape=pop)) +
  xlab(paste("PC1 (",Per_exp[1],"%)",sep="")) + 
  ylab(paste("PC2 (",Per_exp[2],"%)",sep="")) + 
  geom_point(size=3,stroke=1.2) + 
  stat_ellipse(level=0.75,size=1)+
#  scale_shape_manual(name="Pop", labels=unique(pop)[order(unique(pop))], values=shp)+
#  scale_color_manual(name="Pop", labels=unique(pop)[order(unique(pop))], values=cl)+
  scale_shape_manual(name="Pop", labels=unique(pop2)[order(unique(pop2))], values=shp)+
  scale_color_manual(name="Pop", labels=unique(pop2)[order(unique(pop2))], values=cl)+
  labs(color="") + 
  theme_classic()

require(cowplot)
require(ggrepel)
require(seqinr)

