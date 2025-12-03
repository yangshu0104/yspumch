## DEG
library(limma)
library(ggrepel)
library(ggplot2)
library(dplyr)
data<-read.csv('GSE24250_exp.csv',row.names = NULL)
data<-data[!duplicated(data$X),]
rownames(data)<-data$X
data<-data[,-1]
data<-log2(data)
boxplot(data,outline=FALSE,notch=T,las=2)
group<-c(rep('Con',6),rep('HD',8))
group=factor(group,levels = c("Con","HD"))
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(data)
cont.matrix <- makeContrasts(HD - Con , levels = design)
fit<-lmFit(data,design)
fit2<-contrasts.fit(fit,cont.matrix)
fit2<-eBayes(fit2)
tempOutput <- topTable(fit2, coef = 1, n = Inf)
DEG_limma <- na.omit(tempOutput)
logFC = 1
P.Value = 0.05
k1 <- (DEG_limma$P.Value < P.Value) & (DEG_limma$logFC < -logFC)
k2 <- (DEG_limma$P.Value < P.Value) & (DEG_limma$logFC > logFC)
DEG_limma <- mutate(DEG_limma, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_limma$change)
p<-ggplot(data=DEG_limma,aes(x=logFC,y=-log10(P.Value)))+geom_point(alpha=0.4,size=2,aes(color=change))+ylab('-log10(Pvalue)')+scale_color_manual(values = c('blue4','grey','red3'))+geom_vline(xintercept = c(-logFC,logFC),lty=4,col='black',lwd=0.8)+geom_hline(yintercept = -log10(P.Value),lty=4,col='black',lwd=0.8)+theme_bw()
DEG_limma$ID<-rownames(DEG_limma)
up_data<-filter(DEG_limma,change=='up')%>%distinct(ID,.keep_all = TRUE)%>%top_n(10,-log10(P.Value))
down_data<-filter(DEG_limma,change=='down')%>%distinct(ID,.keep_all = TRUE)%>%top_n(10,-log10(P.Value))
p3<-p+geom_point(data=up_data,aes(x=logFC,y=-log10(P.Value)),color='red3',size=4.5,alpha=0.2)+geom_label_repel(data=up_data,aes(x=logFC,y=-log10(P.Value),label=ID),seed = 233,size=3.5,color='black',min.segment.length = 0,force=2,force_pull = 2,box.padding = 0.4,max.overlaps = Inf)+geom_point(data=down_data,aes(x=logFC,y=-log10(P.Value)),color='blue4',size=4.5,alpha=0.2)+geom_label_repel(data=down_data,aes(x=logFC,y=-log10(P.Value),label=ID),seed = 233,size=3.5,color='black',min.segment.length = 0,force=2,force_pull = 2,box.padding = 0.4,max.overlaps = Inf)
data<-read.csv('GSE1751_exp.csv',row.names = NULL)
data<-data[!duplicated(data$X),]
rownames(data)<-data$X
data<-data[,-1]
data<-log2(data)
boxplot(data,outline=FALSE,notch=T,las=2)
group<-c(rep('Con',14),rep('HD',12))
group=factor(group,levels = c("Con","HD"))
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(data)
cont.matrix <- makeContrasts(HD - Con , levels = design)
fit<-lmFit(data,design)
fit2<-contrasts.fit(fit,cont.matrix)
fit2<-eBayes(fit2)
tempOutput <- topTable(fit2, coef = 1, n = Inf)
DEG_limma <- na.omit(tempOutput)
logFC = 1
P.Value = 0.05
k1 <- (DEG_limma$adj.P.Val < P.Value) & (DEG_limma$logFC < -logFC)
k2 <- (DEG_limma$adj.P.Val < P.Value) & (DEG_limma$logFC > logFC)
DEG_limma <- mutate(DEG_limma, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_limma$change)
p<-ggplot(data=DEG_limma,aes(x=logFC,y=-log10(adj.P.Val)))+geom_point(alpha=0.4,size=2,aes(color=change))+ylab('-log10(adj.P.Val)')+scale_color_manual(values = c('blue4','grey','red3'))+geom_vline(xintercept = c(-logFC,logFC),lty=4,col='black',lwd=0.8)+geom_hline(yintercept = -log10(P.Value),lty=4,col='black',lwd=0.8)+theme_bw()
DEG_limma$ID<-rownames(DEG_limma)
up_data<-filter(DEG_limma,change=='up')%>%distinct(ID,.keep_all = TRUE)%>%top_n(10,-log10(adj.P.Val))
down_data<-filter(DEG_limma,change=='down')%>%distinct(ID,.keep_all = TRUE)%>%top_n(10,-log10(adj.P.Val))
p3<-p+geom_point(data=up_data,aes(x=logFC,y=-log10(adj.P.Val)),color='red3',size=4.5,alpha=0.2)+geom_label_repel(data=up_data,aes(x=logFC,y=-log10(adj.P.Val),label=ID),seed = 233,size=3.5,color='black',min.segment.length = 0,force=2,force_pull = 2,box.padding = 0.4,max.overlaps = Inf)+geom_point(data=down_data,aes(x=logFC,y=-log10(adj.P.Val)),color='blue4',size=4.5,alpha=0.2)+geom_label_repel(data=down_data,aes(x=logFC,y=-log10(adj.P.Val),label=ID),seed = 233,size=3.5,color='black',min.segment.length = 0,force=2,force_pull = 2,box.padding = 0.4,max.overlaps = Inf)
p3
