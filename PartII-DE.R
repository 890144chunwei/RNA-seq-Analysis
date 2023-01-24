library("DESeq2")
library(ggplot2)
library(pheatmap)
library(dplyr)


Sample <- c()
Condition <- c()
FC_doxo_all <- data.frame()
for (i in c("NT","Doxo")) {
  for (j in c("wt","het")) {
    for (k in 1:3){
      list <- list()
      name<- paste(i,"_",j,k ,sep="")
      assign(name,read.delim(paste("~/Dropbox/Mac_Desktop/RNAseq_doxo_fc/",i,"_","LSK_",j,k,"_fc.txt", sep = ""), comment.char="#"))
      Sample[length(Sample)+1] <- c((paste(i, paste(j,k,sep = ""),sep="_")))
      Condition[length(Condition)+1] <- c(paste(i,j,sep="_"))
      
    }
  }
}
info_doxo_all <-data.frame(Sample,Condition)


FC_doxo_all <- data.frame(NT_wt1[,7], NT_wt2[,7], NT_wt3[,7], NT_het1[,7], NT_het2[,7], NT_het3[,7], Doxo_wt1[,7], Doxo_wt2[,7], Doxo_wt3[,7], Doxo_het1[,7], Doxo_het2[,7], Doxo_het3[,7])
row.names(FC_doxo_all) <- NT_wt1[,1]



#Prepare DESeq matrix
ddsdoxo_all <- DESeqDataSetFromMatrix(countData = FC_doxo_all, colData = info_doxo_all, design = ~Condition)
ddsdoxo_KO <- DESeqDataSetFromMatrix(countData = FC_doxo_KO, colData = info_doxo_KO, design = ~Condition)
ddsdoxo_NT <- DESeqDataSetFromMatrix(countData = FC_doxo_NT, colData = info_doxo_NT, design = ~Condition)


#remove lowly expressed genes
keepdoxo_all <- rowSums2(counts(ddsdoxo_all)) >= 150
ddsdoxo_all <- ddsdoxo_all[keepdoxo_all,]
keepdoxo_KO <- rowSums2(counts(ddsdoxo_KO)) >= 90
ddsdoxo_KO <- ddsdoxo_KO[keepdoxo_KO,]
keepdoxo_NT <- rowSums2(counts(ddsdoxo_NT)) >= 90
ddsdoxo_NT <- ddsdoxo_NT[keepdoxo_NT,]


#DEseq
DEdoxo_all <- DESeq(ddsdoxo_all)
DEdoxo_allshr <- lfcShrink(DEdoxo_all, coef=2, type='apeglm')
plotDispEsts(DEdoxo_all, ylim =c(1e-4,1e2))
DEdoxo_KO <- DESeq(ddsdoxo_KO)
DEdoxo_KOshr <- lfcShrink(DEdoxo_KO, coef=2, type='apeglm')
plotDispEsts(DEdoxo_KO, ylim =c(1e-3,1e2))
DEdoxo_NT <- DESeq(ddsdoxo_NT)
DEdoxo_NTshr <- lfcShrink(DEdoxo_NT, coef=2, type='apeglm')
plotDispEsts(DEdoxo_NT, ylim =c(1e-4,1e2))


##MA plot for DEseq result
plotMA(results(DEdoxo_all, contrast = c('Condition','NT_het','NT_wt'), alpha = 0.05), ylim =c(-15,15),main="Doxoall: NT_het vs NT_wt")
plotMA(results(DEdoxo_all, contrast = c('Condition','Doxo_wt','NT_wt'), alpha = 0.05), ylim =c(-15,15),main="Doxoall: Doxo_wt vs NT_wt")
plotMA(results(DEdoxo_all, contrast = c('Condition','Doxo_het','NT_het'), alpha = 0.05), ylim =c(-15,15),main="Doxoall: Doxo_het vs NT_het")
plotMA(results(DEdoxo_all, contrast = c('Condition','Doxo_het','Doxo_wt'), alpha = 0.05), ylim =c(-15,15),main="Doxoall: Doxo_het vs Doxo_wt")

#export normalized read counts
normCounts_doxo_all <- counts(DEdoxo_all,normalized = T)
write.csv(normCounts_doxo_all,"~/Desktop/220624_rnaseq_Doxoall_cnt.csv")
Cnt_doxo_all <- read.csv("~/Desktop/220624_rnaseq_Doxoall_cnt.csv", row.names=1)

#PCA analysis
tnormCounts_doxo_all <- t(normCounts_doxo_all)
pca_doxo_all <- prcomp(tnormCounts_doxo_all, scale. = TRUE)

#ggplot for PCA
summary(pca_doxo_all)
pca_doxo_all$x
pca_doxo_all <- data.frame(info_doxo_all[,2], pca_doxo_all$x[,1:2])
colnames(pca_doxo_all) <- c('Condition', 'PC1', 'PC2')
ggplot(pca_doxo_all, aes(PC1, PC2, group=Condition)) + geom_point(size=8, aes(shape = Condition, color = Condition))+ 
  scale_shape_manual(values = c(15,16,15,16))+
  scale_color_manual(values = c("red3","hotpink2", "blue4","steelblue3"))+
  labs(x= "PC1", y= "PC2") + 
  theme(axis.text.y= element_text(size =16), axis.text.x= element_text(size =16), axis.title.x.bottom = element_text(size=20),axis.title.y = element_text(size=20))+
  xlim(-150,150) + ylim(-100, 100)

#DESeq results
res_doxo_all_doxo <- results(DEdoxo_all, contrast = c("Condition", "Doxo_het", "Doxo_wt"), alpha = 0.05, cooksCutoff = FALSE)
res_doxo_all_NT <- results(DEdoxo_all, contrast = c("Condition", "NT_het", "NT_wt"), alpha = 0.05, cooksCutoff = FALSE)

#output DESeq results
res_doxo_all_doxo <- res_doxo_all_doxo[order(res_doxo_all_doxo$padj),]
res_doxo_all_NT <- res_doxo_all_NT[order(res_doxo_all_NT$padj),]
res_doxo_KO <- res_doxo_KO[order(res_doxo_KO$padj),]

write.csv(res_doxo_all_doxo, "~/Desktop/220624_rnaseq_doxoall_doxo_pval.csv")
write.csv(res_doxo_all_NT, "~/Desktop/220624_rnaseq_doxoall_NT_pval.csv")
write.csv(res_doxo_KO, "~/Desktop/220624_rnaseq_doxo_KO_pval.csv")


Pval_all_doxo <- read.csv("~/Desktop/220624_rnaseq_doxoall_doxo_pval.csv", row.names = 1)
Pval_all_NT <- read.csv("~/Desktop/220624_rnaseq_doxoall_NT_pval.csv", row.names = 1)
Pval_doxo_KO <- read.csv("~/Desktop/220624_rnaseq_doxo_KO_pval.csv", row.names = 1)

##Reoragnize DESeq result
Pval_all_doxo$sig <- ifelse(Pval_all_doxo$padj <= 0.05, "yes", "no")
Pval_all_NT$sig <- ifelse(Pval_all_NT$padj <= 0.05, "yes", "no")
Pval_doxo_KO$sig <- ifelse(Pval_doxo_KO$padj <= 0.05, "yes", "no")
Pval_doxo_NT$sig <- ifelse(Pval_doxo_NT$padj <= 0.05, "yes", "no")

##Filter outliers
Pval_all_doxo<- na.omit(Pval_all_doxo)
Pval_all_NT<- na.omit(Pval_all_NT)
Pval_doxo_KO<- na.omit(Pval_doxo_KO)

#remove fold change outlier
Pval_all_doxo <- filter(Pval_all_doxo,log2FoldChange <10, -log2FoldChange<10)
Pval_all_NT <- filter(Pval_all_NT,log2FoldChange <10, -log2FoldChange<10)
Pval_doxo_KO <- filter(Pval_doxo_KO,log2FoldChange <10, -log2FoldChange<10)
Pval_doxo_NT <- filter(Pval_doxo_NT,log2FoldChange <10, -log2FoldChange<10)

#Esembl id to gene name
row.names(Ensembl.id.list) <- Ensembl.id.list[,2]

Pval_all_doxo <- merge(Pval_all_doxo, Ensembl.id.list ,by=0)
row.names(Pval_all_doxo) <- Pval_all_doxo[,9]
Pval_all_doxo <- Pval_all_doxo[,-9:-10]
Pval_all_doxo <- Pval_all_doxo[,-1]
Pval_all_NT <- merge(Pval_all_NT, Ensembl.id.list ,by=0)
row.names(Pval_all_NT) <- Pval_all_NT[,9]
Pval_all_NT <- Pval_all_NT[,-9:-10]
Pval_all_NT <- Pval_all_NT[,-1]
Pval_doxo_KO <- merge(Pval_doxo_KO, Ensembl.id.list ,by=0)
row.names(Pval_doxo_KO) <- Pval_doxo_KO[,9]
Pval_doxo_KO <- Pval_doxo_KO[,-9:-10]
Pval_doxo_KO <- Pval_doxo_KO[,-1]
Pval_doxo_NT <- merge(Pval_doxo_NT, Ensembl.id.list ,by=0)
row.names(Pval_doxo_NT) <- Pval_doxo_NT[,9]
Pval_doxo_NT <- Pval_doxo_NT[,-9:-10]
Pval_doxo_NT <- Pval_doxo_NT[,-1]


write.csv(res_doxo_all_doxo, "~/Desktop/220624_rnaseq_doxoall_doxo_pval.csv")
write.csv(res_doxo_all_NT, "~/Desktop/220624_rnaseq_doxoall_NT_pval.csv")

write.csv(Pval_all_doxo, "~/Desktop/Pval_all_doxo.csv")
write.csv(Pval_all_NT, "~/Desktop/Pval_all_NT.csv")
write.csv(Pval_doxo_KO, "~/Desktop/Pval_doxo_KO.csv")
write.csv(Pval_doxo_NT, "~/Desktop/Pval_doxo_NT.csv")


Pval_all_doxo <- read.csv("~/Desktop/Pval_all_doxo.csv", row.names=1)
Pval_all_NT <- read.csv("~/Desktop/Pval_all_NT.csv", row.names=1)
Pval_doxo_KO <- read.csv("~/Desktop/Pval_doxo_KO.csv", row.names=1)
Pval_doxo_NT <- read.csv("~/Desktop/Pval_doxo_NT.csv", row.names=1)


Pval_all_doxo <- Pval_all_doxo[order(Pval_all_doxo$padj),]
Pval_all_NT <- Pval_all_NT[order(Pval_all_NT$padj),]


##MA plot
ggplot(Pval_nt_het,aes(x= log10(baseMean), y= -log2FoldChange, color = sig)) + 
  geom_point()+ylim(-15, 15) + xlim(1,6.2)+ ggtitle("Doxo: HET vs WT") +
  theme(legend.key.size = (unit(1,'cm')), axis.text.y= element_text(size =16), 
        axis.text.x= element_text(size =16), axis.title.x.bottom = element_text(size=20),
        axis.title.y = element_text(size=20), plot.title= element_text(size=28))

#Volcano plot
ggplot(Pval_all_doxo, aes(x= -log2FoldChange, y= -log10(padj), color = sig)) + geom_point() + 
  ylim(0,10)+ xlim(-5,5) + ggtitle("All_doxo: Mut vs WT") +
  theme(legend.key.size = (unit(1,'cm')), axis.text.y= element_text(size =16), 
        axis.text.x= element_text(size =16), axis.title.x.bottom = element_text(size=20),
        axis.title.y = element_text(size=20), plot.title= element_text(size=28))+
  scale_color_manual(values = c("darkgrey","blue1"))

library(EnhancedVolcano)
EnhancedVolcano(Pval_all_doxo, lab = NA, x='log2FoldChange',y='padj',
                title = 'All_doxo: Mut vs WT', pCutoff = 5e-2, pointSize = 1.5,
                col=c('darkgrey','darkgrey', 'darkgrey','red'), colAlpha = 0.5,
                FCcutoff = 1.2, legendPosition = 'right') + xlim(-10,10)+ ylim(0,15)
EnhancedVolcano(Pval_all_NT, lab = NA, x='log2FoldChange',y='padj',
                title = 'All_NT: Mut vs WT', pCutoff = 5e-2, pointSize = 1.5,
                col=c('darkgrey','darkgrey', 'darkgrey','red'), colAlpha = 0.5,
                FCcutoff = 1.25, legendPosition = 'right') + xlim(-10,10)+ ylim(0,10)

#pheatmap for doxorubicin signature
signi_wt <- subset(Pval_nt_wt, padj <=0.05)
signi_wt <- merge(Cnt_doxo, signi_wt,by=0)
signi_wt_merge <- data.frame(signi_wt[,2:4],signi_wt[,8:10], signi_wt[,5:7],signi_wt[,11:13])
row.names(signi_wt_merge) <- signi_wt$Row.names

signi_wt_merge <- merge(Ensembl.id.list, signi_wt_merge, by=0)
row.names(signi_wt_merge) <- signi_wt_merge$Gene.name
signi_wt_merge <- signi_wt_merge[,5:16]

pheatmap(log2(signi_wt_merge + 1), cluster_cols = F, cellwidth = 24, cellheight = 5.5,border_color = NA ,cutree_rows = 2,treeheight_row = 15, fontsize_col = 0.5,fontsize_row = 5, legend_breaks = c(-3,0,3), scale = 'row')

#pheatmap for doxorubicin het vs wt
signi_doxo_0.05 <- subset(Pval_doxo, padj <=0.05)
signi_doxo_0.05 <- merge(Cnt_doxo, signi_doxo_0.05,by=0)
signi_doxo_merge <- data.frame(signi_doxo_0.05[,8:10], signi_doxo_0.05[,11:13])
row.names(signi_doxo_merge) <- signi_doxo_0.05$Row.names

signi_doxo_0.001_merge <- merge(Ensembl.id.list, signi_doxo_0.001_merge, by=0)
row.names(signi_doxo_0.001_merge) <- signi_doxo_0.001_merge$Gene.name
signi_doxo_0.001_merge <- signi_doxo_0.001_merge[,5:10]

pheatmap(log2(signi_doxo_0.001_merge + 1), cluster_cols = F, cellwidth = 25, cellheight = 6,border_color = NA ,cutree_rows = 2,treeheight_row = 15, fontsize_col = 0.5,fontsize_row = 5, legend_breaks = c(-1.5,0,1.5), scale = 'row')


#GSEA and MsigDB
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(enrichplot)
library(clusterProfiler)
library(ggplot2)

BiocManager::install("org.Mm.eg.db", character.only = TRUE)
BiocManager::install("org.Hs.eg.db", character.only = TRUE)
library("org.Mm.eg.db", character.only = TRUE)
library("org.Hs.eg.db", character.only = TRUE)

Pval_all_doxo_f <- subset(Pval_all_doxo, baseMean >=150)
Genelist_all_doxo <- Pval_all_doxo_f$log2FoldChange
names(Genelist_all_doxo) <- row.names(Pval_all_doxo_f)
Genelist_all_doxo <- na.omit(Genelist_all_doxo)
Genelist_all_doxo <- sort(Genelist_all_doxo, decreasing = TRUE)
Gse_all_doxo <- gseGO(geneList = Genelist_all_doxo, ont = "BP", keyType = "ENSEMBL", 
                 minGSSize = 80, maxGSSize = 400, pvalueCutoff = 1, verbose = TRUE,
                  OrgDb = "org.Mm.eg.db",pAdjustMethod = "none",eps=0)
Gsea_all_doxo <- as.data.frame(Gse_all_doxo)
write.csv(Gsea_all_doxo, "~/Desktop/Gsea_all_doxo.csv")
require(DOSE)
dotplot(Gse_all_doxo, font.size=8, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("GSEA: Doxorubicin Treatment")

Pval_doxo_KO_f <- subset(Pval_doxo_KO, baseMean >=150)
Genelist_doxo_KO <- Pval_doxo_KO_f$log2FoldChange
names(Genelist_doxo_KO) <- row.names(Pval_doxo_KO_f)
Genelist_doxo_KO <- na.omit(Genelist_doxo_KO)
Genelist_doxo_KO <- sort(Genelist_doxo_KO, decreasing = TRUE)
Gse_doxo_KO <- gseGO(geneList = Genelist_doxo_KO, ont = "ALL", keyType = "ENSEMBL", 
                      minGSSize = 80, maxGSSize = 400, pvalueCutoff = 0.05, verbose = TRUE,
                      OrgDb = "org.Mm.eg.db",pAdjustMethod = "BH", eps=0)
Gsea_doxo_KO <- as.data.frame(Gse_doxo_KO)
require(DOSE)
dotplot(Gse_doxo_KO, font.size=8, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("GSEA: Doxorubicin Treatment")

Pval_all_NT_f <- subset(Pval_all_NT, baseMean >=150)
Genelist_all_doxo <- Pval_all_NT_f$log2FoldChange
names(Genelist_all_NT) <- row.names(Pval_all_NT_f)
Genelist_all_NT <- na.omit(Genelist_all_NT)
Genelist_all_NT <- sort(Genelist_all_NT, decreasing = TRUE)
Gse_all_NT <- gseGO(geneList = Genelist_all_NT, ont = "ALL", keyType = "ENSEMBL", 
                      minGSSize = 30, maxGSSize = 300, pvalueCutoff = 0.05, verbose = TRUE,
                      OrgDb = "org.Mm.eg.db",pAdjustMethod = "BH", eps=0)
Gsea_all_NT <- as.data.frame(Gse_all_NT)
require(DOSE)
dotplot(Gse_all_NT, font.size=8, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("GSEA: Doxorubicin Treatment")
rx2wk_B, font.size=8, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("GSEA: Doxorubicin Treatment")

#

install.packages("ggnewscale")
library(ggnewscale)
cnetplot(Gse_trx, categorySize="pvalue", foldChange=Genelist_trx, showCategory = 3)

install.packages("ggridges")
library(ggridges)
ridgeplot(Gse_doxo) + labs(x = "enrichment distribution")+ theme(axis.text.y = element_text(size=6))

#gseaplot(Gse_doxo, by = "all", title = Gse_doxo$Description[1], geneSetID = 1)
gseaplot2(Gse_all_doxo, geneSetID = 211, title = Gse_all_doxo$Description[211])
Gse_all_doxo$NES[211]
Gse_all_doxo$p.adjust[211]
gseaplot2(Gse_all_doxo, geneSetID = 213, title = Gse_all_doxo$Description[213])
Gse_all_doxo$NES[213]
Gse_all_doxo$p.adjust[213]


##Comparing DDR subsets between doxo and transplant
row.names(Cnt_doxo_all) <- gsub(pattern = "[.].*", replacement = "", x= row.names(Cnt_doxo_all))
DE_doxo_DDR <- merge(Genelist_pDNAdamage_m,Cnt_doxo_all, by=0)
row.names(DE_doxo_DDR) <- DE_doxo_DDR[,3]
DE_trx2wk_DDR <- merge(Genelist_pDNAdamage_m,Cnt_trx2wk_LSK, by=0)
row.names(DE_trx2wk_DDR) <- DE_trx2wk_DDR[,3]
DE_doxo_DDR <- DE_doxo_DDR[,-1:-2]
DE_trx2wk_DDR <- DE_trx2wk_DDR[,-1:-2]
DE_doxo_DDR$sumwt <- rowSums(DE_doxo_DDR[,8:9])
DE_doxo_DDR$sumhet <- rowSums(DE_doxo_DDR[,12:13])
DE_doxo_DDR <- filter(DE_doxo_DDR, sumhet - sumwt >1000)
DE_doxo_DDR <- DE_doxo_DDR[,-14:-15]
DE_trx2wk_DDR$sumwt <- rowSums(DE_trx2wk_DDR[,2:3])
DE_trx2wk_DDR$sumhet <- rowSums(DE_trx2wk_DDR[,4:5])
DE_trx2wk_DDR <- filter(DE_trx2wk_DDR, sumhet - sumwt >0)
DE_trx2wk_DDR <- filter(DE_trx2wk_DDR, DE_trx2wk_DDR[,5] - DE_trx2wk_DDR[,3] >0)
DE_trx2wk_DDR <- DE_trx2wk_DDR[,-6:-7]
write.csv(DE_doxo_DDR, "~/Desktop/DE_DDR_doxo.csv")
write.csv(DE_doxo_DDR, "~/Desktop/DE_DDR_trx2wk.csv")

row.names(Cnt_trx) <- gsub(pattern = "[.].*", replacement = "", x= row.names(Cnt_trx))
DE_trx_DDR <- merge(Genelist_pDNAdamage_m,Cnt_trx, by=0)
write.csv(DE_trx_DDR, "~/Desktop/DE_DDR_trx.csv")

DE_DDR_comb <- merge(DE_trx2wk_DDR, DE_doxo_DDR, by=0)
row.names(DE_DDR_comb) <- DE_DDR_comb[,1]
DE_DDR_comb<- DE_DDR_comb[,-1:-2]
which(row.names(DE_DDR_comb) =='Arid1a')
DE_DDR_comb <- DE_DDR_comb[-4,]
pheatmap(log2(DE_DDR_comb[row1,1:4] + 1), cluster_cols = F, cluster_rows = F, cellwidth = 40, cellheight = 14,border_color = NA,treeheight_row = 15, fontsize_col = 5,fontsize_row = 14, legend_breaks = c(-1.5,0,1.5), scale = 'row')
col2<-c(7,13,16,17)
row1<-c(4,12,16,17,11,10,18,1,13,15,5,9,6,7,8,2,3)
pheatmap(log2(DE_DDR_comb[row1,col2] + 1), cluster_cols = F,cluster_rows = F, cellwidth = 40, cellheight = 14,border_color = NA,treeheight_row = 15, fontsize_col = 5,fontsize_row = 14, legend_breaks = c(-1.5,0,1.5), scale = 'row')
col3 <-c(7,6,9,10)
pheatmap(log2(DE_DDR_comb[row1,col3] + 1), cluster_cols = F, cluster_rows = F, cellwidth = 40, cellheight =14,border_color = NA,treeheight_row = 15, fontsize_col = 5,fontsize_row = 14, legend_breaks = c(-1.5,0,1.5), scale = 'row')
write.csv(DE_DDR_comb, "~/Desktop/DE_DDR_comb.csv")

##Comparing stem cell subsets between doxo and transplant
DE_doxo_HSC <- merge(HSC.genelist,Cnt_doxo_all, by=0)
row.names(DE_doxo_HSC) <- DE_doxo_HSC[,2]
DE_trx2wk_HSC <- merge(HSC.genelist,Cnt_trx2wk_LSK, by=0)
row.names(DE_trx2wk_HSC) <- DE_trx2wk_HSC[,2]
DE_doxo_HSC <- DE_doxo_HSC[,-1:-2]
DE_trx2wk_HSC <- DE_trx2wk_HSC[,-1:-2]
DE_doxo_HSC$sumwt <- rowSums(DE_doxo_HSC[,9:10])
DE_doxo_HSC$sumhet <- rowSums(DE_doxo_HSC[,12:13])
DE_doxo_HSC_f <- filter(DE_doxo_HSC, sumhet - sumwt > 0)
DE_doxo_HSC_f <- DE_doxo_HSC_f[,-14:-15]
DE_trx2wk_HSC$sumwt <- rowSums(DE_trx2wk_HSC[,2:3])
DE_trx2wk_HSC$sumhet <- rowSums(DE_trx2wk_HSC[,4:5])
DE_trx2wk_HSC_f <- filter(DE_trx2wk_HSC, sumhet - sumwt >0)
DE_trx2wk_HSC_f <- filter(DE_trx2wk_HSC_f, DE_trx2wk_HSC_f[,5] - DE_trx2wk_HSC_f[,3] >0)
DE_trx2wk_HSC_f <- DE_trx2wk_HSC_f[,-6:-7]
write.csv(DE_doxo_DDR, "~/Desktop/DE_DDR_doxo.csv")
write.csv(DE_HSC_comb, "~/Desktop/DE_HSC_comb.csv")
aaa <- c(13,14,16,17)

DE_HSC_comb <- merge(DEtrx, DE_doxo_HSC, by=0)
row.names(DE_HSC_comb) <- DE_HSC_comb[,1]
DE_HSC_comb<- DE_HSC_comb[,-1:-2]
which(row.names(DE_HSC_comb) =='BCL2L11')
DE_HSC_comb <- DE_HSC_comb[-1,]
pheatmap(log2(DE_HSC_comb[,11:14] + 1), cluster_cols = F, cellwidth = 40, cellheight = 10,border_color = NA,treeheight_row = 15, fontsize_col = 5,fontsize_row = 8, legend_breaks = c(-1.5,0,1.5), scale = 'row')
row2<-c(6,8,1,10,3,12,18,16,17,4,11,7,2,15,5,13,9,14)
pheatmap(log2(DE_HSC_comb[row2,3:6] + 1), cluster_cols = F,cluster_rows = F, cellwidth = 40, cellheight = 14,border_color = NA,treeheight_row = 15, fontsize_col = 2,fontsize_row = 14, legend_breaks = c(-1.5,0,1.5), scale = 'row')
pheatmap(log2(DE_HSC_comb[row2,11:14] + 1), cluster_cols = F, cluster_rows = F, cellwidth = 40, cellheight =14,border_color = NA,treeheight_row = 15, fontsize_col = 5,fontsize_row = 14, legend_breaks = c(-1.5,0,1.5), scale = 'row')
write.csv(DE_DDR_comb, "~/Desktop/DE_DDR_comb.csv")


#DE_DDR_doxo <- DE_DDR[,1:6]
DE_HSC_trx <- DE_HSC[,1:4]
row.names(DE_HSC_trx)<- DE_HSC[,1]
DE_HSC_trx <- filter(DE_HSC_trx, HET1 > WT2, HET1> WT3, HET2> WT2, HET2>WT3)

pheatmap(log2(DE_HSC_trx + 1), cluster_cols = F,cellwidth = 60, cellheight = 10,border_color = NA, fontsize_col = 8,treeheight_row = 15,legend_breaks = c(-1,0,1), fontsize_row = 9, scale = 'row')
##
DE_HSC_doxo <- DE_HSC[,6:11]
DE_HSC_doxo <- DE_HSC_doxo[,-2]
DE_HSC_doxo <- DE_HSC_doxo[,-3]

DE_HSC_doxo$sum_wt <- rowSums(DE_HSC_doxo[,1:2])
DE_HSC_doxo$sum_het <- rowSums(DE_HSC_doxo[,3:4])
DE_HSC_doxo <- filter(DE_HSC_doxo, sum_het - sum_wt>0)
DE_HSC_doxo <- DE_HSC_doxo[,1:4]
DE_HSC_doxo <- filter(DE_HSC_doxo, doxo_het2 > doxo_wt1, doxo_het2> doxo_wt3, doxo_het3> doxo_wt1, doxo_het3>doxo_wt3)

pheatmap(log2(DE_HSC_doxo + 1), cluster_cols = F,cellwidth = 60, cellheight = 10,border_color = NA, fontsize_col = 8,treeheight_row = 15,legend_breaks = c(-1,0,1), fontsize_row = 10, scale = 'row')



#signi_wt <- data.frame(signi_wt[,2:4],signi_wt[,8:10])
#signi_het <- data.frame(signi_wt[,5:7],signi_wt[,11:13])

#row.names(signi_wt) <- signi_wt$Row.names
#row.names(signi_het) <- signi_het$Row.names

#pheatmap(log2(signi_wt + 1), scale = 'row', show_rownames = F)
#pheatmap(log2(signi_het + 1), scale = 'row', show_rownames = F)

## Ensembl id to gene name using biomaRt
#library(biomaRt)
#BiocManager::install("org.Mm.eg.db")
#library(org.Mm.eg.db)

#row.names(Ensembl.id.list) <- Ensembl.id.list[,2]
#row.names(Trx_wt1) <- Trx_wt1[,1]
#test <- merge(Trx_wt1, Ensembl.id.list ,by=0)

