## including manhattan plot and boxplot for the main results

setwd("~/Documents/Projects/Asthma remission/DMA_results/")


####################### complete remission #######################

load("blood/asthma_model2_bloed_complete_remission.txt.Rdata")

anno<-read.csv("~/Documents/Projects/PIAMA/HumanMethylation450.csv",stringsAsFactors = F)
all.anno<-anno[match(rownames(p.all),anno$IlmnID),]

man<-data.frame(SNP=as.character(all.anno$IlmnID),CHR=as.numeric(all.anno$CHR),BP=all.anno$MAPINFO,P=p.all$pvalue)
bonfer<--log10(0.05/dim(man)[1])

intSNP<-"cg24788483"

man.sig<-man[grep(intSNP,man$SNP),]

y<- -log10(man.sig$P) +0.019*12

# calculate the cumulative position of SNP
width<-c()
cum<-c()
cum[1]<-0
for (i in 1:22){
  width[i]<-max(man$BP[which(man$CHR==i)])-min(man$BP[which(man$CHR==i)])
  cum[i+1]<-cum[i]+width[i]
}
x<-c()
for (i in 1:length(y)){
  if(man.sig$CHR[i]<=12)
    x[i]<-man.sig$BP[i]+cum[as.numeric(man.sig$CHR[i])]
  else
    x[i]<-1.5*man.sig$BP[i]+cum[as.numeric(man.sig$CHR[i])]
}

manhattan(man,suggestiveline = F,genomewideline = bonfer,col = c("grey", "gray50"),highlight=intSNP)
annotate.cex=1.5;annotate.font=3
text(x,y,intSNP,cex=annotate.cex,adj=c(0,0.48),font=annotate.font,col="red")

####################### all remission #######################
rm(list=ls())

load("blood/asthma_model2_bloed_all_remission.txt.Rdata")

anno<-read.csv("~/Documents/Projects/PIAMA/HumanMethylation450.csv",stringsAsFactors = F)
all.anno<-anno[match(rownames(p.all),anno$IlmnID),]

man<-data.frame(SNP=as.character(all.anno$IlmnID),CHR=as.numeric(all.anno$CHR),BP=all.anno$MAPINFO,P=p.all$pvalue)
bonfer<--log10(0.05/dim(man)[1])

intSNP<-"cg13378519"
man.sig<-man[grep(intSNP,man$SNP),]

y<- -log10(man.sig$P) +0.019*12

# calculate the cumulative position of SNP
width<-c()
cum<-c()
cum[1]<-0
for (i in 1:22){
  width[i]<-max(man$BP[which(man$CHR==i)])-min(man$BP[which(man$CHR==i)])
  cum[i+1]<-cum[i]+width[i]
}
x<-c()
for (i in 1:length(y)){
  if(man.sig$CHR[i]<=12)
    x[i]<-man.sig$BP[i]+cum[as.numeric(man.sig$CHR[i])]
  else
    x[i]<-1.5*man.sig$BP[i]+cum[as.numeric(man.sig$CHR[i])]
}

manhattan(man,suggestiveline = F,genomewideline = bonfer,col = c("grey", "gray50"),highlight=intSNP)
annotate.cex=1.5;annotate.font=3
text(x,y,intSNP,cex=annotate.cex,adj=c(0,0.48),font=annotate.font,col="red")






