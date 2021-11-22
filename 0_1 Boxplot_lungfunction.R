### -----------------------------------------------------
### boxplot for lung function at baseline and last visit
### for project asthma resmission
### Cancan Qi 31-01-2019
### modified 25-09-2019

library(ggplot2)
library(ggpubr)
rm(list=ls())

setwd("~/Documents/Projects/Asthma remission/")
pheno<-read.csv2("./data/Pheno_baseline_last_all.csv")
pheno$remission.NEW<-as.factor(pheno$remission.NEW)
levels(pheno$remission.NEW)<-c("PerAsthma","CliR_only","ComR")

baseline<-pheno[,c("Samplename","FEV1_B","FVC_B","FEV1pred_B","FVCpred_B","remission.NEW")]
colnames(baseline)<-c("Samplename","FEV1","FVC","FEV1pred","FVCpred","Remission")
baseline$category<-1
last<-pheno[,c("Samplename","FEV1","FVC","FEV1pred","FVCpred","remission.NEW")]
colnames(last)<-c("Samplename","FEV1","FVC","FEV1pred","FVCpred","Remission")
last$category<-2

df<-rbind(baseline,last)
df$category<-as.factor(df$category)
levels(df$category)<-c("Baseline","Last visit")

## FEV1
pdf("./plot/baseline_last_visit/FEV1_star.pdf",width = 8,height = 6)
ggboxplot(df, x = "Remission", y = "FEV1",
               color = "Remission", palette = "jco",
               add = "jitter",
               facet.by = "category", short.panel.labs = FALSE)+
  stat_compare_means(comparisons = list(c("PerAsthma","CliR_only"),
                                        c("PerAsthma","ComR"),
                                        c("CliR_only","ComR")),
                     label = "p.signif",method = "t.test")
dev.off()
# or label = "p.format"

## FVC
pdf("./plot/baseline_last_visit/FVC_star.pdf",width = 8,height = 6)
ggboxplot(df, x = "Remission", y = "FVC",
          color = "Remission", palette = "jco",
          add = "jitter",
          facet.by = "category", short.panel.labs = FALSE)+
  stat_compare_means(comparisons = list(c("PerAsthma","CliR_only"),
                                        c("PerAsthma","ComR"),
                                        c("CliR_only","ComR")),
                     label = "p.signif",method = "t.test")
dev.off()

## FEV1pred
pdf("./plot/baseline_last_visit/FEV1pred_star.pdf",width = 8,height = 6)
ggboxplot(df, x = "Remission", y = "FEV1pred",
          color = "Remission", palette = "jco",
          add = "jitter",
          facet.by = "category", short.panel.labs = FALSE)+
  stat_compare_means(comparisons = list(c("PerAsthma","CliR_only"),
                                        c("PerAsthma","ComR"),
                                        c("CliR_only","ComR")),
                     label = "p.signif",method = "t.test")
dev.off()

## FVCpred
pdf("./plot/baseline_last_visit/FVCpred_star.pdf",width = 8,height = 6)
ggboxplot(df, x = "Remission", y = "FVCpred",
          color = "Remission", palette = "jco",
          add = "jitter",
          facet.by = "category", short.panel.labs = FALSE)+
  stat_compare_means(comparisons = list(c("PerAsthma","CliR_only"),
                                        c("PerAsthma","ComR"),
                                        c("CliR_only","ComR")),
                     label = "p.signif",method = "t.test")
dev.off()
