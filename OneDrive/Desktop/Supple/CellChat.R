## Cellchat
library(CellChat)
library(ComplexHeatmap)
load('srat.Rdata')
Idents(srat) <- 'orig.ident'
HD <- subset(srat, idents = c('GSM8403739_HD_1', 'GSM8403740_HD_2', 'GSM8403741_HD_3'))
WT <- subset(srat, idents = c('GSM8403742_WT_1', 'GSM8403743_WT_2', 'GSM8403744_WT_3'))
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
FeaturePlot(srat,features = c('NRG1'),cols = c('lightgrey','red'),split.by = 'group')
FeaturePlot(srat,features = c('NRG3'),cols = c('lightgrey','red'),split.by = 'group')
FeaturePlot(srat,features = c('ERBB4'),cols = c('lightgrey','red'),split.by = 'group')
VlnPlot(srat,features=c('NRG1'),group.by='group')
VlnPlot(srat,features=c('NRG3'),group.by='group')
VlnPlot(srat,features=c('ERBB4'),group.by='group')

