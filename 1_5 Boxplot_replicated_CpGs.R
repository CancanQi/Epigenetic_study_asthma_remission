## We have two CpGs replicated, and this is to make boxplot of the two CpGs, in both blood and nasal brushes

rm(list=ls())

## load packages
library(ggplot2)

## choose tissue and phenotype
tissue = "bloed"
# tissue = "brush"

## setwd
setwd("~/Documents/Projects/Asthma remission")

samplename<- tolower(colnames(betan))
tissue.index<- grep(tissue,samplename)
beta_matrix<- t(betan[, tissue.index])
rm(betan)
samplename.tissue<- tolower(rownames(beta_matrix))

pheno1<- read.csv("./data/brushforcheng_plus_new.csv")
pheno2<- read.csv("./data/remission_samples_table.csv",sep=",",header=T)

pheno<- cbind(pheno1,pheno2[,c("group.names","age")])
index.tissue<- match(tolower(samplename.tissue),tolower(as.character(pheno$Samplename)))
pheno.tissue<- pheno[index.tissue,c("remission.NEW","age","sex","smoking","packyears","ICS","Samplename")]
colnames(pheno.tissue)<- c("remission","age", "sex","smoking","pack_years","ICS","samplename")
PHENO<- pheno.tissue

beta_matrix <- beta_matrix[(which(tolower(rownames(beta_matrix)) %in% tolower(as.character(PHENO$samplename)))),]
cpgname<- c("cg13378519", "cg24788483")
beta<- data.frame(beta_matrix[,cpgname])
asthma<- PHENO$remission

asthma<- as.character(factor(asthma, labels=c("Asthma","Clin remission","Comp remission")))
beta$asthma_groups<- asthma

g1<- ggplot(beta, aes(x=asthma_groups,y=cg13378519))+ geom_boxplot(aes(fill=asthma_groups)) + ylab("Beta value") +xlab(" ")+ theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) 
g1 + ggtitle("cg13378519 (PEX11B) in nasal brush from discovery cohort") + theme(plot.title = element_text(hjust = 0.5))

g2<- ggplot(beta, aes(x=asthma_groups,y=cg24788483))+ geom_boxplot(aes(fill=asthma_groups)) + ylab("Beta value") +xlab(" ")+ theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) 
g2 + ggtitle("cg24788483 in nasal brush from discovery cohort") + theme(plot.title = element_text(hjust = 0.5))

## in blood, change the tissue at the begining
g3<- ggplot(beta, aes(x=asthma_groups,y=cg13378519))+ geom_boxplot(aes(fill=asthma_groups)) + ylab("Beta value") +xlab(" ")+ theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) 
g3 + ggtitle("cg13378519 (PEX11B) in blood from discovery cohort") + theme(plot.title = element_text(hjust = 0.5))

g4<- ggplot(beta, aes(x=asthma_groups,y=cg24788483))+ geom_boxplot(aes(fill=asthma_groups)) + ylab("Beta value") +xlab(" ")+ theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) 
g4 + ggtitle("cg24788483 in blood from discovery cohort") + theme(plot.title = element_text(hjust = 0.5))


