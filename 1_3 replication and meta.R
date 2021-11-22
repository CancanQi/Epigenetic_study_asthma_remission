## the results of whole blood were repliated in two cohorts, EGEA and Lifelines

#1# ################################# replication in EGEA #################################

# send out: 129 CpGs of complete remission; 7 CpGs of all remission
# model:  glmer( cbind(Meth, Unmeth) ~ MainVar + (1|nofamil) + age + sex + smoking + eosino + neutro + baso + mono + lympho)
# pheno: for 3 phenotype: clinical remission (all remission), complete remission, clinical remission only
# summary the esults of complete remission:

library(stringr)
################## complete remission ##################
setwd("~/Documents/Projects/Asthma remission/Replication/EGEA/15_Dec_EGEA")
comrep<-read.csv("./results/Results_EWAS_pool_Remission_Var_compRem_E2_NoFamily.csv",stringsAsFactors = F,dec=".")
comrep<-comrep[,c("X","Estimate_MainVar","SE_MainVar","Stat_MainVar","Pval_WDT")]

chr<-str_split_fixed(comrep$X,":",2)
comrep$chr<-as.numeric(substring(chr[,1],4,))
comrep$start<-as.numeric(str_split_fixed(chr[,2],"-",2)[,1])
comrep$end<-as.numeric(str_split_fixed(chr[,2],"-",2)[,2])

comsig<-read.table("~/Documents/Projects/Asthma remission/DMA_results/blood/blood_ComRem_genomewide_sig.txt",header = T)

index<-c()
for(i in 1:129){
if(comsig$MAPINFO[i] %in% comrep$start){
index[i]<-which(comsig$CHR[i]==comrep$chr & comsig$MAPINFO[i]==comrep$start)}
else
{index[i]<-NA}
}

df<-cbind(comsig,comrep[index,1:5])
df<-df[order(df$Pval_WDT),]

## 129 CpGs sent for replication, 125 were available in EGEA, 3 CpGs were nominal significant (p<0.05), within the 3, only one 
## with the same direction (cg24788483)

## meta analysis using weighted Z-score
beta1<-df$BETA
beta2<-df$Estimate_MainVar
pval1<-df$P_VAL
pval2<-df$Pval_WDT

n1<-44+10
n2<-109
## or for each probe the number can be different

nsum <- n1+n2
w1 <- sqrt(n1) / sqrt(nsum)
w2 <- sqrt(n2) / sqrt(nsum)

z1<-c()
for (i in 1:(length(pval1)))
{
  pvali<-pval1[i]
  betai<-beta1[i]
  if (is.na(pvali)==FALSE){
if ( betai > 0 ) {
  zi <- qnorm( pvali / 2 );
} else { 
  zi <- -(qnorm( pvali / 2 ));
}}
  else
  {zi<-NA}
  z1<-c(z1,zi)
}

z2<-c()
for (i in 1:(length(pval2)))
{
  pvali<-pval2[i]
  betai<-beta2[i]
  if (is.na(pvali)==FALSE){
    if ( betai > 0 ) {
      zi <- qnorm( pvali / 2 );
    } else { 
      zi <- -(qnorm( pvali / 2 ));
    }}
  else
  {zi<-NA}
  z2<-c(z2,zi)
}

zsum <- ( w1*z1 ) + ( w2*z2 ) 

pmeta <- pnorm(-(abs(zsum))) * 2

################## all remission (clinical remission) ##################
rm(list=ls())
allrep<-read.csv("./results/Results_",stringsAsFactors = F,dec=".")
allrep<-allrep[,c("X","Estimate_MainVar","SE_MainVar","Stat_MainVar","Pval_WDT")]

chr<-str_split_fixed(allrep$X,":",2)
allrep$chr<-as.numeric(substring(chr[,1],4,))
allrep$start<-as.numeric(str_split_fixed(chr[,2],"-",2)[,1])
allrep$end<-as.numeric(str_split_fixed(chr[,2],"-",2)[,2])

allsig<-read.table("~/Documents/Projects/Asthma remission/DMA_results/blood/blood_AllRem_genomewide_sig.txt",header = T)

index<-c()
for(i in 1:129){
if(allsig$MAPINFO[i] %in% allrep$start){
index[i]<-which(allsig$CHR[i]==allrep$chr & allsig$MAPINFO[i]==allrep$start)}
else
{index[i]<-NA}
}

df<-cbind(allsig,allrep[index[1:7],1:5])
df<-df[order(df$Pval_WDT),]

## 7 CpGs sent for replication, 7 were available at EGEA, 1 was nominal significant (p<0.05), none with the same direction


#2# ################################# replication in LifeLines #################################
rm(list=ls())
setwd("~/Documents/Projects/Asthma remission/Replication/lifelines")

allrem<-read.table("clinremVSactive_houseman.txt",header = T)
allsig<-read.table("~/Documents/Projects/Asthma remission/DMA_results/blood/blood_AllRem_genomewide_sig.txt",header = T)

df<-allrem[match(allsig$probeID,allrem$CpG),]
## 7 CpGs sent for replication, 7 were availbale in Lifelines, 2 was nominal significant (p<0.05), 1 with the same direction











    
