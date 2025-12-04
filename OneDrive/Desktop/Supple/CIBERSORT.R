## CIBERSORT
library(reshape2)
library(ggpubr)
library(ggplot2)
library(Hmisc)
library(pheatmap)
library(dplyr)
library(limma)
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
DEG <- DEG_limma[abs(DEG_limma$logFC)>1 & DEG_limma$P.Value<0.05,]
DEG1 <- data[rownames(DEG),]
write.table(DEG1,'DEG.txt',sep = '\t',row.names = T,col.names = T,quote = F)
source('./CIBERSORT.R')
LM22<-'LM22.txt'
DEG<-'DEG.txt'
result1 <- CIBERSORT('LM22.txt','DEG.txt', perm = 1000, QN = T)
celldata<-as.data.frame(result1[,1:22])
celldata$group<-c(rep('Con',6),rep('HD',8))
celldata$sample<-row.names(celldata)
celldata_new<-melt(celldata)
colnames(celldata_new) = c("Group", "Sample", "Celltype", "Composition")
ggplot(celldata_new, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell Proportion", x =  NULL, title = "Cell Proportion")+  
  geom_boxplot(aes(fill = Group), position = position_dodge(0.5), width = 0.5, outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#096EA9", "#B33D27")) +
  theme_bw() + 
  theme(plot.title = element_text(size = 12,color="black",hjust = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)) + 
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)
mygene <- c('ADRA2A')
t<-celldata[,1:22]
t<-t(t)
nc = t(rbind(t,data[mygene,]))
m = rcorr(nc)$r[1:nrow(t),(ncol(nc)-length(mygene)+1):ncol(nc)]
p = rcorr(nc)$P[1:nrow(t),(ncol(nc)-length(mygene)+1):ncol(nc)]
m<-na.omit(m)
p<-na.omit(p)
tmp <- paste(case_when(as.vector(p) < 0.01 ~ "**",
                       as.vector(p) < 0.05 ~ "*",
                       TRUE ~ ""), nrow = nrow(p))
p1 <- pheatmap(t(m),
               display_numbers =t(tmp),
               angle_col =90,
               color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
               border_color = "white",
               cellwidth = 20, 
               cellheight = 20,
               width = 7, 
               height=9.1,
               treeheight_col = 0,
               treeheight_row = 0,cluster_rows = F,cluster_cols = F,fontsize_col = 10, ,fontsize_number = 14)

