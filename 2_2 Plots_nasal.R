

### scater plot for cross tissue effect
### this is to make scatter plot of effect size of CpGs in blood and nasal

library(ggplot2)

setwd("~/Documents/Projects/Asthma remission/")
load("./DMA_results/blood/asthma_model2_bloed_all_remission.txt.Rdata")
brush<-read.table("./sva_brush/be/asthma_model1_brush_all_remission.txt.gz")

load("./DMA_results/blood/asthma_model2_bloed_complete_remission.txt.Rdata")
brush<-read.table("./sva_brush/be/asthma_model1_brush_complete_remission.txt.gz")

sig<-p.all[which(p.all$pvalue<1e-6),]

sig.brush<-brush[match(rownames(sig),brush$probeID),]

df<-data.frame(blood_coef=sig$coef,brush_coef=sig.brush$BETA)

df$direction<-0

df$direction[which(df$blood_coef*df$brush_coef>0)]<-1

df$direction<-as.factor(df$direction)

ggplot(df,aes(x=blood_coef, 
              y=brush_coef,color=direction))+geom_point()

