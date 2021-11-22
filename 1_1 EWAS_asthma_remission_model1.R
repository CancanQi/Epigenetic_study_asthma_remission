## this script is to run EWAS for both complete remission and all remission for both blood and nasal brushes
## without cell type adjustment for blood samples
## the analysis was performed on peregrine
## modified from Cheng's code 
## 2017/4/

## load pacakges
rm(list=ls())

library(compare)
library(parallel)
library(plyr)
library(MASS)
library(sandwich)
library(lmtest)
library(data.table)
library(R.utils)

## pheno and tissue

tissue = "bloed"
remission = "all_remission"

## output file

filename<- paste0("asthma_model1_",tissue,"_",remission,".txt")

## load data
setwd("/home/p282473/projects/asthma_remission")
load("./data/beta_BBRMI.Rdata")
samplename<- tolower(colnames(betan))
tissue.index<- grep(tissue,samplename)
beta_matrix<- t(betan[, tissue.index])
samplename.tissue<- tolower(rownames(beta_matrix))

## read phenotype
pheno1<- read.csv("./Phenotype/brushforcheng_plus_new.csv")
pheno2<- read.csv("./Phenotype/remission_samples_table.csv",sep=",",header=T)

## deal with batch effect
batch <- c()
batch.values <- apply(pheno2,1,function(x) {
  if(substring(x[1],1,1)=="2") {
    batch<- c(batch, 1)
  } else {
    batch<- c(batch, 0)
  }
})

## phenotype restructure

pheno<- cbind(pheno1,pheno2[,c("group.names","age")],batch.values)
index.tissue<- match(tolower(samplename.tissue),tolower(as.character(pheno$Samplename)))
pheno.tissue<- pheno[index.tissue,c("remission.NEW","age","sex","smoking","packyears","batch.values","Samplename")]
colnames(pheno.tissue)<- c("remission","age", "sex","smoking","pack_years","batch","samplename")
PHENO<- pheno.tissue

## select all remission or complete remission
## all remission is complete remission and only clinical remission

if (remission  == "all_remission") {
  PHENO$remission[PHENO$remission == 2] <- 1	
} else {
  PHENO <- PHENO[which(PHENO$remission!=1),]
  rownames(beta_matrix)
  beta_matrix <- beta_matrix[(which(tolower(rownames(beta_matrix)) %in% tolower(as.character(PHENO$samplename)))),]
}

## make sure the sample order is consistent
if(compare(tolower(as.character(PHENO$samplename)),tolower(rownames(beta_matrix)))[[1]]==FALSE) stop("Sample order does not match.")
PHENO<- droplevels(PHENO)

## main function of rubost linear regression model
RLMtest = function(meth_matrix,methcol,asthma, X1, X2, X3, X4,X5) {
  mod = try(rlm(meth_matrix[,methcol]~asthma+X1+X2+X3+X4+X5,maxit=200))
  if(class(mod) == "try-error"){
    print(paste("error thrown by column", methcol))
    invisible(rep(NA, 3))
  }else cf = coeftest(mod, vcov=vcovHC(mod, type="HC0"))
  cf[2, c("Estimate", "Std. Error", "Pr(>|z|)")]
}

system.time(ind.res <- mclapply(mc.preschedule=T,mc.cores=1,setNames(seq_len(ncol(beta_matrix)), dimnames(beta_matrix)[[2]]), 
                                RLMtest, meth_matrix=beta_matrix, asthma=PHENO[,1], X1=PHENO[,2], X2=PHENO[,3], X3=PHENO[,4], X4=as.numeric(as.character(PHENO[,5])), X5=PHENO[,6]))

setattr(ind.res, 'class', 'data.frame')
setattr(ind.res, "row.names", c(NA_integer_,4))
setattr(ind.res, "names", make.names(names(ind.res), unique=TRUE))
probelistnames <- names(ind.res)
all.results <- t(data.table(ind.res))
all.results<-data.table(all.results)
all.results[, probeID := probelistnames]
setnames(all.results, c("BETA","SE", "P_VAL", "probeID")) # rename columns
setcolorder(all.results, c("probeID","BETA","SE", "P_VAL"))
rm(probelistnames, ind.res)

# export table of results

write.table(all.results, filename,na="NA")
gzip(filename)


