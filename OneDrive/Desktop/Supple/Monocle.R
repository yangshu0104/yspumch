## Monocle
library(monocle)
load('srat_hd.Rdata')
scRNA = subset(srat,downsample = 2000)
ct <- scRNA@assays$RNA$counts
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
expressed_genes<-VariableFeatures(scRNA)
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
library(ggsci)
p1 = plot_cell_trajectory(sc_cds)+ scale_color_nejm()
p2 = plot_cell_trajectory(sc_cds, color_by = 'Pseudotime') 
p3 = plot_cell_trajectory(sc_cds, color_by = 'celltype')  + scale_color_npg()
BEAM_res<-BEAM(sc_cds[deg$gene_short_name,], branch_point = 1, cores = 2, progenitor_method = 'duplicate')
BEAM_res<-BEAM_res[order(BEAM_res$qval),]
BEAM_res<-BEAM_res[,c('gene_short_name','pval','qval')]
p <- plot_genes_branched_heatmap(sc_cds[row.names(subset(BEAM_res,qval<1e-4)),],branch_point = 1,num_clusters = 3,cores = 1,show_rownames = F,return_heatmap = T)
write.csv(p$annotation_row,file = 'clustering_hd.csv')
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
markers <- read.csv(file = "clustering_hd.csv", stringsAsFactors = FALSE,row.names = 1)
markers$ID<-rownames(markers)
gid <- bitr(unique(markers$ID), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
markers <- full_join(markers, gid, by=c('ID' = 'SYMBOL'))
markers<-na.omit(markers)
markers_2<-markers[markers$Cluster=='1' | markers$Cluster=='2',]
GO = compareCluster(ENTREZID ~ Cluster, data = markers_2, fun='enrichGO',OrgDb='org.Hs.eg.db')
dotplot(GO, label_format=40) + theme(axis.text.x = element_text(angle=45, hjust=1)) + scale_color_gradient(high="#4b5cc4",low="#FE8D3C")
load('srat_wt.Rdata')
scRNA = subset(srat,downsample = 2000)
ct <- scRNA@assays$RNA$counts
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
expressed_genes<-VariableFeatures(scRNA)
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
p1 = plot_cell_trajectory(sc_cds)+ scale_color_nejm()
p2 = plot_cell_trajectory(sc_cds, color_by = 'Pseudotime') 
p3 = plot_cell_trajectory(sc_cds, color_by = 'celltype')  + scale_color_npg()
BEAM_res<-BEAM(sc_cds[deg$gene_short_name,], branch_point = 1, cores = 2, progenitor_method = 'duplicate')
BEAM_res<-BEAM_res[order(BEAM_res$qval),]
BEAM_res<-BEAM_res[,c('gene_short_name','pval','qval')]
p <- plot_genes_branched_heatmap(sc_cds[row.names(subset(BEAM_res,qval<1e-4)),],branch_point = 1,num_clusters = 3,cores = 1,show_rownames = F,return_heatmap = T)
write.csv(p$annotation_row,file = 'clustering_wt.csv')
data <- read.csv(file = "clustering_wt.csv", stringsAsFactors = FALSE,row.names = 1)
df<-bitr(unique(rownames(data)),fromType = "SYMBOL",toType = c('ENTREZID'),OrgDb = org.Hs.eg.db)
data$symbol<-rownames(data)
DEG<-merge(data,df,by.y="SYMBOL",by.x='symbol')
DEG1<-DEG[DEG$Cluster=='1',]
DEG1=DEG1[,'ENTREZID']
enrich.go<-enrichGO(gene = DEG1,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID',ont = 'ALL',pvalueCutoff = 1,readable = T)
barplot(enrich.go)
DEG2<-DEG[DEG$Cluster=='2',]
DEG2=DEG2[,'ENTREZID']
enrich.go<-enrichGO(gene = DEG2,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID',ont = 'ALL',pvalueCutoff = 1,readable = T)
barplot(enrich.go)

