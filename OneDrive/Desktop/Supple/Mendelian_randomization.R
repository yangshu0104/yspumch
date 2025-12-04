## Mendelian randomization
library(TwoSampleMR)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
data1 <- fread('finngen_R12_G6_NEUATR.gz')
data1<-rename(data1,SNP=rsids)
data <- fread('17682_1_CD46_CD46.txt.gz')
outTab<-subset(data, Pval<5e-08)
write.csv(outTab, file="exposure.pvalue.csv", row.names=F)
exposure_dat<-read_exposure_data(filename = 'exposure.pvalue.csv',sep = ",",snp_col = "rsids",beta_col = "Beta",se_col = "SE",effect_allele_col = "effectAllele",other_allele_col = "otherAllele",pval_col = "Pval", eaf_col = "ImpMAF", samplesize_col = "N",clump = F)
exposure_dat$exposure<-'CD46'
exposure_dat_clumped <- clump_data(exposure_dat, clump_kb=10000, clump_r2=0.01)
intersection_dat <- merge(exposure_dat_clumped,data1,by.x="SNP",by.y="SNP")
write.csv(intersection_dat,file="intersection_dat.csv")
outcome_dat <- read_outcome_data(snps = exposure_dat$SNP,filename = "intersection_dat.csv", sep = ",",snp_col = "SNP",beta_col = "beta",se_col = "sebeta",effect_allele_col = "alt",other_allele_col = "ref",pval_col = "pval",eaf_col = "af_alt")
outcome_dat$outcome<-'HD'
dat<-harmonise_data(exposure_dat=exposure_dat_clumped,outcome_dat=outcome_dat)
mrResult=mr(dat)
mrTab=generate_odds_ratios(mrResult)
heterTab=mr_heterogeneity(dat)
pleioTab=mr_pleiotropy_test(dat)
res_single=mr_singlesnp(dat)
mr_forest_plot(res_single)
mr_scatter_plot(mr_results=mrResult, dat=dat)
mr_funnel_plot(singlesnp_results = res_single)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
data<-read.csv('MR.csv',row.names = 1)
data$Sig <- case_when(as.numeric(data$pval) < 0.0001~"****",
                      as.numeric(data$pval) < 0.001~"***",
                      as.numeric(data$pval) < 0.01~"**",
                      as.numeric(data$pval) < 0.05~"*",
                      T~ "")
range(data$or_lci95)
range(data$or_uci95)
data$Group <- factor(rownames(data),levels = rev(rownames(data)))
p1 <-  ggplot(data,aes(x = or, y = Group)) +
  geom_point(color = "white") + 
  labs(title = "MR Analysis", x = "Odd ratio",y = NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,2),
                     breaks = c(0.5,1,1.5),labels = c(0.5,1,1.5)) + 
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 8-0.5, ymax = 8+0.5, fill = '#f3eded', alpha = .9) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 7-0.5, ymax = 7+0.5, fill = '#eaeded', alpha = .9) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 6-0.5, ymax = 6+0.5, fill = '#f3eded', alpha = .9) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 5-0.5, ymax = 5+0.5, fill = '#eaeded', alpha = .9) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 4-0.5, ymax = 4+0.5, fill = '#f3eded', alpha = .9) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 3-0.5, ymax = 3+0.5, fill = '#eaeded', alpha = .9) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 2-0.5, ymax = 2+0.5, fill = '#f3eded', alpha = .9) +
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 1-0.5, ymax = 1+0.5, fill = '#e1e4e4', alpha = .9)
mytheme <- theme_bw() +
  theme(plot.title = element_text(size = 16,hjust = 0.5,face = "bold",color = "#951515"),
        axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16,face = "bold",color = "#951515"), 
        axis.title.y = element_blank(), 
        axis.ticks.length = unit(0.7,"lines"), 
        axis.ticks = element_line(size = 2),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.position = "bottom",
        legend.background = element_rect(fill = 'transparent'))
p3 <- p1 + mytheme + 
  geom_errorbar(aes(xmin = or_lci95 , xmax = or_uci95), 
                position = position_dodge(0.8), 
                width = 0.1, size = 1.5,
                color = "black") +
  geom_point(aes(color = Group,fill = Group),
             shape = 22,size = 10) +
  geom_vline(aes(xintercept = 1),linetype = 'dashed',
             color = 'grey', size = 1.5) +
  scale_color_manual(values = c("#ff6666","#ffc857","#93a8ac","#23949b",'#E41A1C','#802168','#4DAF4A','#377EB8')) +
  scale_fill_manual(values = c("#ff6666","#ffc857","#93a8ac","#23949b",'#E41A1C','#802168','#4DAF4A','#377EB8'))
p3 + geom_text(aes(label = Sig),vjust = 1.5,size = 7)
library(shapviz)
library(xgboost)
library(pROC)
data =read.csv("shap_GSE24250.csv",header = T,check.names = F)
model_xgboost = xgboost(
  data = as.matrix(data[,-1]),
  label = data$Trait,
  max_depth = 3, 
  eta = 1, 
  nthread = 2, 
  nrounds = 10,
  objective = "binary:logistic")
shap_xgboost <- shapviz(model_xgboost,X_pred =data.matrix(data[,-1]))
sv_importance(shap_xgboost,fill = '#0085FF')+theme_bw()
data$pred <- predict(model_xgboost, as.matrix(data[,-1]))
ran_roc<-roc(data$Trait,data$pred,ci=TRUE)
plot(ran_roc,grid = c(0.5,0.2), print.auc = T, print.auc.pattern = 'AUC=%.3f', col = "#ff6666", grid.col = "gray", xlab = "1 - Specificity", ylab = "Sensitivity", font.lab = 2,  cex.axis = 1.2,print.auc.x = 0.6, print.auc.y = 0.8)
data =read.csv("shap_GSE1751.csv",header = T,check.names = F)
model_xgboost = xgboost(
  data = as.matrix(data[,-1]),
  label = data$Trait,
  max_depth = 3, 
  eta = 1, 
  nthread = 2, 
  nrounds = 10,
  objective = "binary:logistic")
data$pred <- predict(model_xgboost, as.matrix(data[,-1]))
ran_roc1<-roc(data$Trait,data$pred,ci=TRUE)
plot(ran_roc1,add=TRUE,col='#23949b',print.auc = T, print.auc.pattern = 'AUC=%.3f',print.auc.x = 0.6, print.auc.y = 0.7)
legend('bottomright',legend = c('GSE24250','GSE1751'),col = c("#ff6666",'#23949b'),lty = 1)

