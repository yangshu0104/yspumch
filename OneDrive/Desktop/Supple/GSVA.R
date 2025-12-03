## GSVA
library(msigdbr)
library(GSVA)
library(dplyr)
library(limma)
library(pheatmap)
GO_df = msigdbr(species = "Homo sapiens",category = "C5") %>% 
  dplyr::select(gene_symbol,gs_exact_source,gs_subcat,gs_name)
GO_df = GO_df[GO_df$gs_subcat!="HPO",]              
GO_df <- GO_df[, c("gs_name", "gene_symbol")]
go_list<-split(GO_df$gene_symbol,GO_df$gs_name)
exprSet <- as.matrix(data) 
gsva_mat <- gsvaParam(exprSet,go_list)
gsva<-gsva(gsva_mat)
group<-c(rep('Con',14),rep('HD',12))
group<-factor(group,levels = c('Con','HD'))
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group)) 
rownames(design) = colnames(data)
cont.matrix <- makeContrasts(HD-Con, levels = design)
fit<-lmFit(gsva,design)
fit2<-contrasts.fit(fit,cont.matrix)
fit2<-eBayes(fit2)
tempOutput <- topTable(fit2, coef = 1, n = Inf)
degs <- na.omit(tempOutput)
padj_cutoff=0.05
log2FC_cutoff=log2(2)
keep <- rownames(degs[degs$adj.P.Val < padj_cutoff & abs(degs$logFC)>log2(2), ])
length(keep)
dat <- gsva[keep[1:25],]
annotation_col<-data.frame(group=group)
rownames(annotation_col)<-colnames(data)
pheatmap::pheatmap(dat, 
                   fontsize_row = 10,
                   height = 10,
                   width=18,
                   annotation_col = annotation_col,
                   show_colnames = F,
                   show_rownames = T,
                   filename = paste0('GSVA_go_heatmap.pdf'))

