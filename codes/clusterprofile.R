#Import your dataset
#Col1 (Genes ID in Entrez form)
#Col2 (LogFoldChange)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("DOSE")
BiocManager::install("enrichplot")


#Access libraries
library(DOSE)
library(enrichplot)
library(ReactomePA)
library(clusterProfiler)

#Run Enrichment
nomes <- geneList1$ID
valores <- geneList1$LogFC
names(valores) <- nomes
genevector <- names(valores)

edo <- enrichPathway(genevector)
egoCC <- enrichGO(genevector, ont="CC",  OrgDb = "org.Hs.eg.db")
egoBP <- enrichGO(genevector, ont="BP",  OrgDb = "org.Hs.eg.db")
egoMF <- enrichGO(genevector, ont="MF",  OrgDb = "org.Hs.eg.db")

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
egoxCC <- setReadable(egoCC, 'org.Hs.eg.db', 'ENTREZID')
egoxBP <- setReadable(egoBP, 'org.Hs.eg.db', 'ENTREZID')
egoxMF <- setReadable(egoMF, 'org.Hs.eg.db', 'ENTREZID')


###############################################################
## PLOTS (Pahtways)

#Dot plot of enriched terms.
dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")

#Network plot of enriched terms.
p1 <- cnetplot(edox, foldChange=valores)
p2 <- cnetplot(edox, foldChange=valores, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p2, ncol=1, rel_widths=c(1, 1))

#Labelling nodes by selected subset
p3 <- cnetplot(edox, node_label="none")
p4 <- cnetplot(edox, node_label="gene") 
p5 <- cnetplot(edox, node_label="category") 
p6 <- cnetplot(edox, node_label="all") 
cowplot::plot_grid(p6, ncol=1)

#Heatmap plot of enriched terms
p7 <- heatplot(edox)
p8 <- heatplot(edox, foldChange=valores)
cowplot::plot_grid(p8, ncol=1)

#Plot for results obtained from hypergeometric test and gene set enrichment analysis.
p9 <- emapplot(edo)
p10 <- emapplot(edo, pie_scale=1.5)
p11 <- emapplot(edo,layout="kk")
p12 <- emapplot(edo, pie_scale=1.5,layout="kk") 
cowplot::plot_grid(p12, ncol=1)



###########################################################
## PLOTS (CC)

#Dot plot of enriched terms.
dotplot(egoCC, showCategory=30) + ggtitle("dotplot for CC")

#Network plot of enriched terms.
p1 <- cnetplot(egoxCC, foldChange=valores)
p2 <- cnetplot(egoxCC, foldChange=valores, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p2, ncol=1, rel_widths=c(1, 1))

#Labelling nodes by selected subset
p3 <- cnetplot(egoxCC, node_label="none")
p4 <- cnetplot(egoxCC, node_label="gene") 
p5 <- cnetplot(egoxCC, node_label="category") 
p6 <- cnetplot(egoxCC, node_label="all") 
cowplot::plot_grid(p6, ncol=1)

#Heatmap plot of enriched terms
p7 <- heatplot(egoxCC)
p8 <- heatplot(egoxCC, foldChange=valores)
cowplot::plot_grid(p8, ncol=1)

#Plot for results obtained from hypergeometric test and gene set enrichment analysis.
p9 <- emapplot(egoCC)
p10 <- emapplot(egoCC, pie_scale=1.5)
p11 <- emapplot(egoCC,layout="kk")
p12 <- emapplot(egoCC, pie_scale=1.5,layout="kk") 
cowplot::plot_grid(p12, ncol=1)


###############################################################
## PLOTS (BP)

#Dot plot of enriched terms.
dotplot(egoBP, showCategory=30) + ggtitle("dotplot for BP")

#Network plot of enriched terms.
p1 <- cnetplot(egoxBP, foldChange=valores)
p2 <- cnetplot(egoxBP, foldChange=valores, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p2, ncol=1, rel_widths=c(1, 1))

#Labelling nodes by selected subset
p3 <- cnetplot(egoxBP, node_label="none")
p4 <- cnetplot(egoxBP, node_label="gene") 
p5 <- cnetplot(egoxBP, node_label="category") 
p6 <- cnetplot(egoxBP, node_label="all") 
cowplot::plot_grid(p6, ncol=1)

#Heatmap plot of enriched terms
p7 <- heatplot(egoxBP)
p8 <- heatplot(egoxBP, foldChange=valores)
cowplot::plot_grid(p8, ncol=1)

#Plot for results obtained from hypergeometric test and gene set enrichment analysis.
p9 <- emapplot(egoBP)
p10 <- emapplot(egoBP, pie_scale=1.5)
p11 <- emapplot(egoBP,layout="kk")
p12 <- emapplot(egoBP, pie_scale=1.5,layout="kk") 
cowplot::plot_grid(p12, ncol=1)


###############################################################
## PLOTS (MF)

#Dot plot of enriched terms.
dotplot(egoMF, showCategory=30) + ggtitle("dotplot for MF")

#Network plot of enriched terms.
p1 <- cnetplot(egoxMF, foldChange=valores)
p2 <- cnetplot(egoxMF, foldChange=valores, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p2, ncol=1, rel_widths=c(1, 1))

#Labelling nodes by selected subset
p3 <- cnetplot(egoxMF, node_label="none")
p4 <- cnetplot(egoxMF, node_label="gene") 
p5 <- cnetplot(egoxMF, node_label="category") 
p6 <- cnetplot(egoxMF, node_label="all") 
cowplot::plot_grid(p6, ncol=1)

#Heatmap plot of enriched terms
p7 <- heatplot(egoxMF)
p8 <- heatplot(egoxMF, foldChange=valores)
cowplot::plot_grid(p8, ncol=1)

#Plot for results obtained from hypergeometric test and gene set enrichment analysis.
p9 <- emapplot(egoMF)
p10 <- emapplot(egoMF, pie_scale=1.5)
p11 <- emapplot(egoMF,layout="kk")
p12 <- emapplot(egoMF, pie_scale=1.5,layout="kk") 
cowplot::plot_grid(p12, ncol=1)


#######################################################################
#Export tables

install.packages("writexl")
library(writexl)

tabelaORA <- edox@result
tabelaCC <- egoxCC@result
tabelaBP <- egoxBP@result
tabelaMF <- egoxMF@result


write_xlsx(tabelaORA, "C:\\Users/otavi/tabelaORA.xlsx")
write_xlsx(tabelaCC, "C:\\Users/otavi/tabelaCC.xlsx")
write_xlsx(tabelaBP, "C:\\Users/otavi/tabelaBP.xlsx")
write_xlsx(tabelaMF, "C:\\Users/otavi/tabelaMF.xlsx")
