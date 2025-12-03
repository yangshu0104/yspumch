## GSEA
library(clusterProfiler)
library(org.Hs.eg.db)
genelist <- DEG_limma$logFC
names(genelist) <- rownames(DEG_limma)
genelist <- na.omit(genelist)
genelist <- sort(genelist,decreasing = TRUE)
gse <- gseGO(geneList = genelist,ont = "ALL",OrgDb = org.Hs.eg.db,keyType = "SYMBOL",pvalueCutoff = 0.05,pAdjustMethod = "none",verbose = TRUE)
require(DOSE)
dotplot(gse,showCategory=5,split=".sign",font.size=11)+facet_grid(.~.sign)
