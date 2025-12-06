## SCENIC
library(SCopeLoomR)
library(AUCell) 
library(SCENIC)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(optparse)
subcell.id<-sample(colnames(srat),800)
scRNAsub<-srat[,subcell.id]
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
closeAllConnections()
gc()
runGenie3(exprMat_filtered_log, scenicOptions, nParts = 5)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1
scenicOptions@settings$seed <- 123
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget"))
library(doParallel)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
regulonAUC <-loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "celltype"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "group"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
regulonActivity_byGroup <- sapply(split(rownames(cellInfo), cellInfo$group),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup), center = T, scale=T))
pheatmap(regulonActivity_byCellType, scale = 'row',show_colnames=T,cluster_rows = F)

