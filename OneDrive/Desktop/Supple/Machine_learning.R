## Machine learning
library(glmnet)
library(dplyr)
data<-read.csv('GSE24250_GSE1751_overlap.csv',row.names = 1)
data<-scale(data)
data<-t(data)
data<-as.data.frame(data)
data$Trait<-'Unknown'
data$Trait[rownames(data)=='normal_1' | rownames(data)=='normal_2' | rownames(data)=='normal_3' | rownames(data)=='normal_4' | rownames(data)=='normal_5' | rownames(data)=='normal_6']<-'0'
data$Trait[rownames(data)=='HD_1' | rownames(data)=='HD_2' | rownames(data)=='HD_3' | rownames(data)=='HD_4' | rownames(data)=='HD_5' | rownames(data)=='HD_6' | rownames(data)=='HD_7' | rownames(data)=='HD_8']<-'1'
data<-select(data,Trait,everything())
data<-as.matrix(data)
mod_cv<-cv.glmnet(data[,-1],data[,1],family='binomial',alpha=1)
plot(mod_cv)
best_lambda<-mod_cv$lambda.min
best_model<-glmnet(data[,-1],data[,1],family = 'binomial',alpha = 1,lambda = best_lambda)
coef(best_model)
write.csv(as.matrix(coef(best_model)),file='lasso.csv')
library(randomForest)
data<-read.csv('GSE24250_GSE1751_overlap.csv',row.names = 1)
data<-t(data)
data<-as.data.frame(data)
data$Trait<-'Unknown'
data$Trait[rownames(data)=='normal_1' | rownames(data)=='normal_2' | rownames(data)=='normal_3' | rownames(data)=='normal_4' | rownames(data)=='normal_5' | rownames(data)=='normal_6']<-'0'
data$Trait[rownames(data)=='HD_1' | rownames(data)=='HD_2' | rownames(data)=='HD_3' | rownames(data)=='HD_4' | rownames(data)=='HD_5' | rownames(data)=='HD_6' | rownames(data)=='HD_7' | rownames(data)=='HD_8']<-'1'
data<-select(data,Trait,everything())
data$Trait<-as.numeric(data$Trait)
colnames(data) <- gsub("-", "_", colnames(data))
set.seed(123456789)
train_randomforest<-randomForest(Trait~.,data=data,ntree=500,importance=TRUE,proximity=TRUE)
plot(train_randomforest)
varImpPlot(train_randomforest,n.var=10)
result<-rfcv(data[,-1],data[,1],cv.fold = 5,scale = 'log',step = 0.99)
with(result,plot(n.var,error.cv,log='x',type='o',lwd=2))
library(e1071)
data<-read.csv('GSE24250_GSE1751_overlap.csv',row.names = 1)
data<-scale(data)
data<-t(data)
data<-as.data.frame(data)
data$Trait<-'Unknown'
data$Trait[rownames(data)=='normal_1' | rownames(data)=='normal_2' | rownames(data)=='normal_3' | rownames(data)=='normal_4' | rownames(data)=='normal_5' | rownames(data)=='normal_6']<-'0'
data$Trait[rownames(data)=='HD_1' | rownames(data)=='HD_2' | rownames(data)=='HD_3' | rownames(data)=='HD_4' | rownames(data)=='HD_5' | rownames(data)=='HD_6' | rownames(data)=='HD_7' | rownames(data)=='HD_8']<-'1'
data<-select(data,Trait,everything())
data$Trait<-as.numeric(data$Trait)
source("./msvmRFE.R")
nfold = 10
nrows = nrow(data)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, data, k=10, halve.above=100)
top.features = WriteFeatures(results, data, save=F)
featsweep = lapply(1:5, FeatSweep.wrap, results, data)
no.info = min(prop.table(table(data[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
PlotErrors(errors, no.info=no.info)
library(VennDiagram)
venn <- read.csv('GSE24250_venn_genes.csv')
a<-venn$Lasso[1:4]
b<-venn$RandomForest[1:3]
c<-venn$SVM[1:2]
venn_list<-list(Lasso=a,randomForest=b,SVM=c)
venn.diagram(venn_list, cex=1.5,fill=c('#FFFFCC','#CCFFFF','#FFCCCC'),cat.col='black',cat.cex=1.5,cat.fontfamily='serif',col=c('#FFFFCC','#CCFFFF','#FFCCCC'),fontfamily='serif',filename = "test2.png",cat.dist=c(0.055,0.055,0.055),cat.pos=c(-27,27,-180))
library(ggpubr)
library(ggplot2)
GSE24250<-read.csv('GSE24250_ADRA2A.csv')
GSE1751<-read.csv('GSE1751_ADRA2A.csv')
GSE24250$Type<-'GSE24250'
GSE1751$Type<-'GSE1751'
GSE<-rbind(GSE24250,GSE1751)
ggplot2::ggplot(GSE,aes(x=Type,y=ADRA2A,fill=Group))+geom_boxplot(alpha=0.7)+scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+theme(panel.grid = element_blank())+stat_compare_means(aes(group = Group),label = "p.signif",method = "wilcox.test",hide.ns = T)+theme_classic()
library(pROC)
GSE24250_ADRA2A<-roc(GSE24250$Group,GSE24250$ADRA2A,levels=c('Con','HD'))
GSE1751_ADRA2A<-roc(GSE1751$Group,GSE1751$ADRA2A,levels=c('Con','HD'))
plot(GSE24250_ADRA2A,grid = c(0.5,0.2), print.auc = T, print.auc.pattern = 'AUC=%.3f', col = "#DD7123", grid.col = "gray", xlab = "False Positive Rate (1 - Specificity)", ylab = "True Positive Rate (Sensitivity)", font.lab = 2,  cex.axis = 1.2,print.auc.x = 0.6, print.auc.y = 0.8)
plot(GSE1751_ADRA2A,add=TRUE,col='blue',print.auc = T, print.auc.pattern = 'AUC=%.3f',print.auc.x = 0.6, print.auc.y = 0.7)
legend('bottomright',legend = c('GSE24250','GSE1751'),col = c("#DD7123",'blue'),lty = 1)

