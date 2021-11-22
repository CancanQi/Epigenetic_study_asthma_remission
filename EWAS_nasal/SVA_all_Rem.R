###########################################
#### SVA to adjust for hidden factors
#### this script is used to generate svas
#### you need add svas to your linear model
#### Cancan Qi 13/02/2018
###########################################

rm(list=ls())

library(sva)
library(MASS)
library(compare)

setwd("/home/p282473/projects/asthma_remission/")
load("data/beta_BBRMI.Rdata")
tissue = "brush"
remission="all_remission"
filename<- paste0("sva_",tissue,"_",remission,".txt")
samplename<- tolower(colnames(betan))
tissue.index<- grep(tissue,samplename)
beta_matrix<- t(betan[, tissue.index])
rm(betan)
samplename.tissue<- tolower(rownames(beta_matrix))

pheno1<- read.csv("./Phenotype/brushforcheng_plus_new(nospace).csv")
pheno2<- read.csv("./Phenotype/remission_samples_table.csv",sep=",",header=T)

batch <- c()
batch.values <- apply(pheno2,1,function(x) {
  if(substring(x[1],1,1)=="2") {
    batch<- c(batch, 1)
  } else {
    batch<- c(batch, 0)
  }
})

pheno<- cbind(pheno1,pheno2[,c("group.names","age")],batch.values)
index.tissue<- match(tolower(samplename.tissue),tolower(as.character(pheno$Samplename)))
pheno.tissue<- pheno[index.tissue,c("remission.NEW","age","sex","smoking","packyears","ICS","batch.values","Samplename")]
colnames(pheno.tissue)<- c("remission","age", "sex","smoking","pack_years","ICS","batch","samplename")
rownames(pheno.tissue)<-pheno.tissue$samplename
PHENO<- pheno.tissue[,-8]

if (remission  == "all_remission") {
  PHENO$remission[PHENO$remission == 2] <- 1	
} else {
  PHENO <- PHENO[which(PHENO$remission!=1),]
  rownames(beta_matrix)
  beta_matrix <- beta_matrix[(which(tolower(rownames(beta_matrix)) %in% tolower(as.character(PHENO$samplename)))),]
}

PHENO<- droplevels(PHENO)
PHENO<-na.omit(PHENO)
beta_matrix<-beta_matrix[match(rownames(PHENO),rownames(beta_matrix)),]
if(compare(tolower(rownames(PHENO)),tolower(rownames(beta_matrix)))[[1]]==FALSE) stop("Sample order does not match.")

beta<-t(beta_matrix)

mod1<-model.matrix(~.,PHENO)
mod0<-model.matrix(~.,PHENO[,-1])
n.sv <- num.sv(dat =beta, mod = mod1, method = "be", vfilter = 2000, B = 1000, seed =1)  ## method of "be"
# n.sv <- num.sv(dat = beta, mod = mod1, method = "leek")   ## method of "leek"
svobj1= sva(beta,mod1,mod0,n.sv=n.sv)
SVs = as.data.frame(svobj1$sv)
colnames(SVs) <-paste0("sv",1:ncol(SVs))
modSv1 = cbind(PHENO,SVs)

setwd("./sva")
save(modSv1,file = "./sva_all_remission.Rdata")



