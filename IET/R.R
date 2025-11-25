## DEG
library(AnnoProbe)
library(GEOquery)
library(stringr)
library(limma)
gset=AnnoProbe::geoChina('GSE24250')
gset=gset[[1]]
exprSet=exprs(gset)
phenoDat=pData(gset)
group_list=ifelse(str_detect(phenoDat$title,"Huntington's"),"HD","control")
group_list=factor(group_list,levels = c("control","HD"))
exprSet <- log2(exprSet)
(gpl=gset@annotation)
ids<-idmap(gpl ,type = 'soft')
exprSet<-filterEM(exprSet,ids)
exprSet$symbol=rownames(exprSet)
exprSet=data.frame(exprSet[-grep("/",exprSet$symbol),])
exprSet=exprSet[,1:14]
design <- model.matrix(~0 + factor(group_list))
colnames(design) = levels(factor(group_list))
rownames(design) = colnames(exprSet)
cont.matrix <- makeContrasts(HD-control, levels = design)
fit<-lmFit(exprSet,design)
fit2<-contrasts.fit(fit,cont.matrix)
fit2<-eBayes(fit2)
tempOutput <- topTable(fit2, coef = 1, n = Inf)
DEG_limma <- na.omit(tempOutput)
DEG_limma <- subset(DEG_limma,abs(logFC)>1 | P.Value < 0.05)

## GSEA
library(msigdbr)
library(dplyr)
library(clusterProfiler)
GO_df = msigdbr(species = "Homo sapiens",category = "C5") %>% 
  dplyr::select(gene_symbol,gs_exact_source,gs_subcat,gs_name)
GO_df = GO_df[GO_df$gs_subcat!="HPO",]              
ge = DEG_limma$logFC
names(ge) = rownames(DEG_limma)
ge = sort(ge,decreasing = T)
em <- GSEA(ge, TERM2GENE = GO_df[, c("gs_name", "gene_symbol")],pvalueCutoff = 1,pAdjustMethod = 'BH')

## GSVA
library(GSVA)
GO_df <- GO_df[, c("gs_name", "gene_symbol")]
go_list<-split(GO_df$gene_symbol,GO_df$gs_name)
exprSet <- as.matrix(exprSet) 
gsva_mat <- gsvaParam(exprSet,go_list)
gsva<-gsva(gsva_mat)

## Machine learning
library(glmnet)
mod_cv<-cv.glmnet(data[,-1],data[,1],family='binomial',alpha=1)
plot(mod_cv)
best_lambda<-mod_cv$lambda.min
best_model<-glmnet(data[,-1],data[,1],family = 'binomial',alpha = 1,lambda = best_lambda)
coef(best_model)
write.csv(as.matrix(coef(best_model)),file='lasso.csv')
library(randomForest)
set.seed(123456789)
train_randomforest<-randomForest(Trait~.,data=data,ntree=500,importance=TRUE,proximity=TRUE)
plot(train_randomforest)
varImpPlot(train_randomforest,n.var=20)
result<-rfcv(data[,-1],data[,1],cv.fold = 5,scale = 'log',step = 0.99)
with(result,plot(n.var,error.cv,log='x',type='o',lwd=2))
library(e1071)
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
library(pROC)
ROC<-roc(data$Group,data$ADRA2A,levels=c('normal','HD'))

## CIBERSORT
library(reshape2)
library(ggpubr)
library(ggplot2)
library(Hmisc)
library(pheatmap)
write.table(exprSet,'data.txt',sep = '\t',row.names = T,col.names = T,quote = F)
source('./CIBERSORT.R')
LM22<-'LM22.txt'
data<-'data.txt'
result1 <- CIBERSORT('LM22.txt','data.txt', perm = 1000, QN = T)
celldata<-as.data.frame(result1[,1:22])
celldata$group<-c('Con','HD','HD','HD','HD','HD','HD','Con','Con','Con','HD','Con','HD','Con')
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
nc = t(rbind(t,exprSet[mygene,]))
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

## Mendelian randomization
library(TwoSampleMR)
library(data.table)
library(tidyr)
data1 <- fread('finngen_R12_G6_NEUATR.gz')
data1<-rename(data1,SNP=rsids)
data <- fread('13405_61_SPINK2_ISK2.txt.gz')
outTab<-subset(data, Pval<5e-08)
write.csv(outTab, file="exposure.pvalue.csv", row.names=F)
exposure_dat<-read_exposure_data(filename = 'exposure.pvalue.csv',sep = ",",snp_col = "rsids",beta_col = "Beta",se_col = "SE",effect_allele_col = "effectAllele",other_allele_col = "otherAllele",pval_col = "Pval", eaf_col = "ImpMAF", samplesize_col = "N",clump = F)
exposure_dat$exposure<-'SPINK2'
exposure_dat_clumped <- clump_data(exposure_dat, clump_kb=10000, clump_r2=0.01)
intersection_dat <- merge(exposure_dat_clumped,data1,by.x="SNP",by.y="SNP")
write.csv(intersection_dat,file="intersection_dat.csv")
outcome_dat <- read_outcome_data(snps = exposure_dat$SNP,filename = "intersection_dat.csv", sep = ",",snp_col = "SNP",beta_col = "Beta",se_col = "SE",effect_allele_col = "Test_Allele",other_allele_col = "Other_Allele",pval_col = "P.val",eaf_col = "Test_Allele_Frequency")
outcome_dat$outcome<-'HD'
dat<-harmonise_data(exposure_dat=exposure_dat_clumped,outcome_dat=outcome_dat)
mrResult=mr(dat)
mrTab=generate_odds_ratios(mrResult)
library(shapviz)
library(xgboost)
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

## scRNA
library(Seurat)
library(harmony)
library(monocle)
samples<-list.files('./GSE271852')
seurat_list<-list()
for (sample in samples) {
  data.path <- paste0("./GSE271852/", sample)
  seurat_data <- Read10X(data.dir = data.path)
  seurat_obj <- CreateSeuratObject(counts = seurat_data,project = sample,min.features = 200,
                                   min.cells = 3)
  seurat_list <- append(seurat_list, seurat_obj)
}
seurat_combined <- merge(seurat_list[[1]], 
                         y = seurat_list[-1],
                         add.cell.ids = samples)
srat<-JoinLayers(seurat_combined)
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
VlnPlot(srat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
srat <- subset(srat, subset = nCount_RNA < 20000 & nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
srat <- NormalizeData(srat)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
srat <- ScaleData(srat)
srat <- RunPCA(srat, features = VariableFeatures(object = srat))
srat <- RunHarmony(srat,group.by.vars='orig.ident',plot_covergence=TRUE)
srat <- RunUMAP(srat, dims = 1:20, reduction = 'harmony')
srat <- FindNeighbors(srat, dims = 1:20,reduction = 'harmony')
srat <- FindClusters(srat, resolution = 0.5)
new.cluster.ids <- c( "Microglial cell", "Microglial cell","Microglial cell", "Oligodendrocyte", "Astrocyte", "Neuron","Oligodendrocyte","Oligodendrocyte","Microglial cell", "Oligodendrocyte", "Astrocyte")
names(new.cluster.ids) <- levels(srat)
srat <- RenameIdents(srat, new.cluster.ids)
DimPlot(srat, label = TRUE, pt.size = 0.5,repel = TRUE) + NoLegend()
DEG<-FindMarkers(srat,ident.1 = 'Disease',ident.2 = 'Normal',group.by = 'group',test.use='wilcox')

## SCENIC
library(SCENIC)
library(AUCell)
library(BiocParallel)
library(doParallel)
library(foreach)
exprMat<-as.matrix(scRNAsub@assays$RNA$data)
cellInfo=data.frame(scRNAsub@meta.data)
colnames(cellInfo)[which(colnames(cellInfo)=="orig.ident")]="sample"
colnames(cellInfo)[which(colnames(cellInfo)=="seurat_clusters")]="cluster"
cellInfo=cellInfo[,c("sample","cluster","celltype","group")]
saveRDS(cellInfo,file="./int/cellinfo.rds")
org='hgnc'
data("defaultDbNames")
dbs <- defaultDbNames[[org]]
myDatasetTitle='HD'
dbDir<-'cisTarget_databases'
motifAnnotations_hgnc <- motifAnnotations
scenicOptions <- initializeScenic(org = org,dbDir = dbDir,dbs = dbs,datasetTitle = myDatasetTitle,nCores = 10)
scenicOptions@inputDatasetInfo$cellinfo <- "./int/cellinfo.rds"
saveRDS(scenicOptions,file = 'int/scenicOptions.Rds')
genesKept <- geneFiltering(exprMat, scenicOptions,minCountsPerGene = 3 * 0.01 * ncol(exprMat),  minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 20)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget"))
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
regulonAUC <-loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
celltype = as.data.frame(subset(cellInfo,select = 'celltype'))
pheatmap(regulonActivity_byCellType, scale = 'row',show_colnames=F,annotation_col = celltype,cluster_rows = F)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled[1:15,], name="Regulon activity")
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "Group", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "celltype"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)

## Monocle
library(monocle)
srat$celltype = Idents(srat)
ct <- srat@assays$RNA$counts
gene_ann <- data.frame(
  gene_short_name = row.names(ct), 
  row.names = row.names(ct)
)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
pd <- new("AnnotatedDataFrame",
          data=scRNA@meta.data)
sc_cds <- newCellDataSet(
  ct, 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds <- estimateSizeFactors(sc_cds)
sc_cds <- estimateDispersions(sc_cds)
expressed_genes<-VariableFeatures(srat)
diff_test_res <- differentialGeneTest(sc_cds[expressed_genes,],
                                      fullModelFormulaStr = " ~ celltype + orig.ident", 
                                      reducedModelFormulaStr = " ~ orig.ident", 
                                      cores = 4)
deg <- subset(diff_test_res,qval<0.01)
deg<-deg[order(deg$qval,decreasing = F),]
deg<-deg[1:1000,]
ordergene<-deg$gene_short_name
sc_cds <- setOrderingFilter(sc_cds, ordergene)
plot_ordering_genes(sc_cds)
sc_cds <- reduceDimension(sc_cds,residualModelFormulaStr = "~orig.ident")
sc_cds <- orderCells(sc_cds)
p1 = plot_cell_trajectory(sc_cds, color_by = 'Pseudotime') 
p2 = plot_cell_trajectory(sc_cds, color_by = 'celltype')
p3 = plot_cell_trajectory(sc_cds, color_by = 'State')
BEAM_res<-BEAM(sc_cds[deg$gene_short_name,], branch_point = 1, cores = 2, progenitor_method = 'duplicate')
BEAM_res<-BEAM_res[order(BEAM_res$qval),]
BEAM_res<-BEAM_res[,c('gene_short_name','pval','qval')]
p <- plot_genes_branched_heatmap(sc_cds[row.names(subset(BEAM_res,qval<1e-4)),],branch_point = 1,num_clusters = 3,cores = 1,show_rownames = F,return_heatmap = T)

## Cellchat
library(CellChat)
library(ComplexHeatmap)
Idents(srat) <- 'orig.ident'
WT <- subset(srat, idents = c('GSM8403742', 'GSM8403743', 'GSM8403744'))
HD <- subset(srat, idents = c('GSM8403739', 'GSM8403740', 'GSM8403741'))
cco.wt <- createCellChat(WT@assays$RNA$data, meta = WT@meta.data, group.by = "celltype")
cco.hd <- createCellChat(HD@assays$RNA$data, meta = HD@meta.data, group.by = 'celltype')
cco.wt <- setIdent(cco.wt, ident.use = 'celltype')
cco.hd <- setIdent(cco.hd, ident.use = 'celltype')
cellchat <- cco.wt
cellchat@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cco.wt <- cellchat
cellchat <- cco.hd
cellchat@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cco.hd <- cellchat
cco.list <- list(wt=cco.wt, hd=cco.hd)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
pathways.show <- c("NRG")
weight.max <- getMaxWeight(cco.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(cco.list)[i]))
}
p <- netVisual_bubble(cellchat, sources.use = c(1:4), targets.use = c(2,3),  comparison = c(1, 2), angle.x = 45)


