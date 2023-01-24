library("DESeq2")
library(ggplot2)
library(pheatmap)
library(dplyr)

Doxo_het1 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_doxo_fc/Doxo_LSK_het1_fc.txt", comment.char="#")
Doxo_het2 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_doxo_fc/Doxo_LSK_het2_fc.txt", comment.char="#")
Doxo_het3 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_doxo_fc/Doxo_LSK_het3_fc.txt", comment.char="#")
Doxo_wt1 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_doxo_fc/Doxo_LSK_wt1_fc.txt", comment.char="#")
Doxo_wt2 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_doxo_fc/Doxo_LSK_wt2_fc.txt", comment.char="#")
Doxo_wt3 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_doxo_fc/Doxo_LSK_wt3_fc.txt", comment.char="#")
NT_wt1 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_doxo_fc/NT_LSK_wt1_fc.txt", comment.char="#")
NT_wt2 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_doxo_fc/NT_LSK_wt2_fc.txt", comment.char="#")
NT_wt3 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_doxo_fc/NT_LSK_wt3_fc.txt", comment.char="#")
NT_het1 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_doxo_fc/NT_LSK_het1_fc.txt", comment.char="#")
NT_het2 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_doxo_fc/NT_LSK_het2_fc.txt", comment.char="#")
NT_het3 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_doxo_fc/NT_LSK_het3_fc.txt", comment.char="#")

Trx2wk_LSK_wt1 <- read.delim("~/Dropbox/Mac_Desktop/Trx2wk_FC/LSK_WT1_fc.txt", comment.char="#")
Trx2wk_LSK_wt3 <- read.delim("~/Dropbox/Mac_Desktop/Trx2wk_FC/LSK_WT3_fc.txt", comment.char="#")
Trx2wk_LSK_mut1 <- read.delim("~/Dropbox/Mac_Desktop/Trx2wk_FC/LSK_Mut1_fc.txt", comment.char="#")
Trx2wk_LSK_mut2 <- read.delim("~/Dropbox/Mac_Desktop/Trx2wk_FC/LSK_Mut2_fc.txt", comment.char="#")
Trx2wk_B_wt1 <- read.delim("~/Dropbox/Mac_Desktop/Trx2wk_FC/Bcell_WT1_fc.txt", comment.char="#")
Trx2wk_B_wt2 <- read.delim("~/Dropbox/Mac_Desktop/Trx2wk_FC/Bcell_WT2_fc.txt", comment.char="#")
Trx2wk_B_mut1 <- read.delim("~/Dropbox/Mac_Desktop/Trx2wk_FC/Bcell_Mut1_fc.txt", comment.char="#")
Trx2wk_B_mut3 <- read.delim("~/Dropbox/Mac_Desktop/Trx2wk_FC/Bcell_Mut3_fc.txt", comment.char="#")
Trx2wk_T_wt1 <- read.delim("~/Dropbox/Mac_Desktop/Trx2wk_FC/Tcell_WT1_fc.txt", comment.char="#")
Trx2wk_T_wt2 <- read.delim("~/Dropbox/Mac_Desktop/Trx2wk_FC/Tcell_WT2_fc.txt", comment.char="#")
Trx2wk_T_wt3 <- read.delim("~/Dropbox/Mac_Desktop/Trx2wk_FC/Tcell_WT3_fc.txt", comment.char="#")
Trx2wk_T_mut1 <- read.delim("~/Dropbox/Mac_Desktop/Trx2wk_FC/Tcell_Mut1_fc.txt", comment.char="#")
Trx2wk_T_mut3 <- read.delim("~/Dropbox/Mac_Desktop/Trx2wk_FC/Tcell_Mut3_fc.txt", comment.char="#")

Trx_wt1 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_trx_fc/Trx_LSK_wt1_fc.txt", comment.char="#")
Trx_wt2 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_trx_fc/Trx_LSK_wt2_fc.txt", comment.char="#")
Trx_wt3 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_trx_fc/Trx_LSK_wt3_fc.txt", comment.char="#")
Trx_het1 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_trx_fc/Trx_LSK_het1_fc.txt", comment.char="#")
Trx_het2 <- read.delim("~/Dropbox/Mac_Desktop/RNAseq_trx_fc/Trx_LSK_het2_fc.txt", comment.char="#")

info_doxo_all <- data.frame(X = c("NT_wt1","NT_wt2","NT_wt3", "NT_het1","NT_het2","NT_het3",
                                  "Doxo_wt1","Doxo_wt2","Doxo_wt3","Doxo_het1","Doxo_het2","Doxo_het3"), 
                          Condition = c("NT_wt","NT_wt","NT_wt","NT_het","NT_het","NT_het",
                                        "Doxo_wt","Doxo_wt","Doxo_wt","Doxo_het","Doxo_het","Doxo_het") )
info_doxo_KO <- data.frame(X = c("Doxo_wt1","Doxo_wt2","Doxo_wt3","Doxo_het1","Doxo_het2","Doxo_het3"), 
                            Condition = c("WT","WT","WT","Mut","Mut","Mut") )
info_doxo_NT <- data.frame(X = c("NT_wt1","NT_wt2","NT_wt3", "NT_het1","NT_het2","NT_het3"), 
                            Condition = c("WT","WT","WT","Mut","Mut","Mut"))
info_trx2wk_LSK <- data.frame(X = c("Trx2wk_LSK_wt1","Trx2wk_LSK_wt3","Trx2wk_LSK_mut1","Trx2wk_LSK_mut2"), 
                           Condition = c("WT","WT","Mut","Mut") )
info_trx2wk_B <- data.frame(X = c("Trx2wk_B_wt1","Trx2wk_B_wt2","Trx2wk_B_mut1","Trx2wk_B_mut3"), 
                              Condition = c("WT","WT","Mut","Mut") )
info_trx2wk_T <- data.frame(X = c("Trx2wk_T_wt1","Trx2wk_T_wt2","Trx2wk_T_wt3","Trx2wk_T_mut1","Trx2wk_T_mut3"), 
                            Condition = c("WT","WT","WT","Mut","Mut") )
info_trx <- data.frame(X = c("Trx_wt1","Trx_wt2","Trx_wt3", "Trx_het1","Trx_het2"), 
                            Condition = c("WT","WT","WT","Mut","Mut") )
info_trx2wk_LSK <- data.frame(X = c("NT_wt1","NT_wt2","NT_wt3", "NT_het1","NT_het2","NT_het3", "Trx2wk_LSK_wt1","Trx2wk_LSK_wt3","Trx_wt1", "Trx2wk_LSK_mut1","Trx2wk_LSK_mut2","Trx_het1"), 
                              Condition = c("WT","WT","WT","Mut","Mut","Mut","WT_Trx","WT_Trx","WT_Trx","Mut_Trx","Mut_Trx","Mut_Trx") )

info_trx2wk_all <- data.frame(X = c("Trx2wk_LSK_wt1","Trx2wk_LSK_wt3","Trx2wk_LSK_mut2","Trx2wk_LSK_mut1","Trx2wk_B_wt1","Trx2wk_B_wt2","Trx2wk_B_mut1","Trx2wk_B_mut3","Trx2wk_T_wt1","Trx2wk_T_wt2","Trx2wk_T_wt3","Trx2wk_T_mut1","Trx2wk_T_mut3"), 
                              Condition = c("WT_LSK","WT_LSK","Mut_LSK","Mut_LSK","WT_B","WT_B","Mut_B","Mut_B","WT_T","WT_T","WT_T","Mut_T","Mut_T" ) )
info_trx_all <- data.frame(X = c("Trx2wk_LSK_wt1","Trx2wk_LSK_wt3","Trx2wk_LSK_mut2","Trx2wk_LSK_mut1","Trx_wt1","Trx_wt2","Trx_wt3", "Trx_het1","Trx_het2"), 
                           Condition = c("WT_2wk","WT_2wk","Mut_2wk","Mut_2wk","WT_5mo","WT_5mo","WT_5mo","Mut_5mo","Mut_5mo") )
info_fig4 <- data.frame(X = c("NT_wt1","NT_wt2","NT_het2","NT_het3","Doxo_wt1","Doxo_wt3","Doxo_het2","Doxo_het3","Trx_wt2","Trx_wt3","Trx_het1","Trx_het2"),
                            Condition = c("NT_wt","NT_wt","NT_het","NT_het","Doxo_wt","Doxo_wt","Doxo_het","Doxo_het","Trx_wt","Trx_wt","Trx_het","Trx_het") )


FC_doxo_all <- data.frame(NT_wt1[,7], NT_wt2[,7], NT_wt3[,7], NT_het1[,7], NT_het2[,7], NT_het3[,7], Doxo_wt1[,7], Doxo_wt2[,7], Doxo_wt3[,7], Doxo_het1[,7], Doxo_het2[,7], Doxo_het3[,7])
row.names(FC_doxo_all) <- NT_wt1[,1]
FC_doxo_KO <- data.frame(Doxo_wt1[,7], Doxo_wt2[,7], Doxo_wt3[,7], Doxo_het1[,7], Doxo_het2[,7], Doxo_het3[,7])
row.names(FC_doxo_KO) <- Doxo_wt1[,1]
FC_doxo_NT <- data.frame(NT_wt1[,7], NT_wt2[,7], NT_wt3[,7], NT_het1[,7], NT_het2[,7], NT_het3[,7])
row.names(FC_doxo_NT) <- NT_wt1[,1]
FC_trx2wk_LSK <- data.frame(Trx2wk_LSK_wt1[,7], Trx2wk_LSK_wt3[,7], Trx2wk_LSK_mut1[,7], Trx2wk_LSK_mut2[,7])
row.names(FC_trx2wk_LSK) <- Trx2wk_LSK_wt1[,1]
FC_trx2wk_B <- data.frame(Trx2wk_B_wt1[,7], Trx2wk_B_wt2[,7], Trx2wk_B_mut1[,7], Trx2wk_B_mut3[,7])
row.names(FC_trx2wk_B) <- Trx2wk_B_wt1[,1]
FC_trx2wk_T <- data.frame(Trx2wk_T_wt1[,7], Trx2wk_T_wt2[,7],Trx2wk_T_wt3[,7], Trx2wk_T_mut1[,7], Trx2wk_T_mut3[,7])
row.names(FC_trx2wk_T) <- Trx2wk_T_wt1[,1]
FC_trx <- data.frame(Trx_wt1[,7], Trx_wt2[,7], Trx_wt3[,7],Trx_het1[,7], Trx_het2[,7])
row.names(FC_trx) <- Trx_wt1[,1]
FC_trx2wk_all <- data.frame(Trx2wk_LSK_wt1[,7], Trx2wk_LSK_wt3[,7], Trx2wk_LSK_mut2[,7], Trx2wk_LSK_mut1[,7],Trx2wk_B_wt1[,7], Trx2wk_B_wt2[,7], Trx2wk_B_mut1[,7], Trx2wk_B_mut3[,7],Trx2wk_T_wt1[,7], Trx2wk_T_wt2[,7],Trx2wk_T_wt3[,7], Trx2wk_T_mut1[,7], Trx2wk_T_mut3[,7])
row.names(FC_trx2wk_all) <- Trx2wk_LSK_wt1[,1]
FC_trx_all <- data.frame(Trx2wk_LSK_wt1[,7], Trx2wk_LSK_wt3[,7], Trx2wk_LSK_mut2[,7], Trx2wk_LSK_mut1[,7],Trx_wt1[,7], Trx_wt2[,7], Trx_wt3[,7],Trx_het1[,7], Trx_het2[,7])
row.names(FC_trx_all) <- Trx2wk_LSK_wt1[,1]
FC_fig_4 <- data.frame(NT_wt1[,7], NT_wt2[,7], NT_het2[,7], NT_het3[,7], Doxo_wt1[,7], Doxo_wt3[,7], Doxo_het2[,7], Doxo_het3[,7], Trx_wt2[,7], Trx_wt3[,7],Trx_het1[,7], Trx_het2[,7])
row.names(FC_fig_4) <- NT_wt1[,1]


#Prepare DESeq matrix
ddsdoxo_all <- DESeqDataSetFromMatrix(countData = FC_doxo_all, colData = info_doxo_all, design = ~Condition)
ddsdoxo_KO <- DESeqDataSetFromMatrix(countData = FC_doxo_KO, colData = info_doxo_KO, design = ~Condition)
ddsdoxo_NT <- DESeqDataSetFromMatrix(countData = FC_doxo_NT, colData = info_doxo_NT, design = ~Condition)
ddstrx2wk_LSK <- DESeqDataSetFromMatrix(countData = FC_trx2wk_LSK, colData = info_trx2wk_LSK, design = ~Condition)
ddstrx2wk_B <- DESeqDataSetFromMatrix(countData = FC_trx2wk_B, colData = info_trx2wk_B, design = ~Condition)
ddstrx2wk_T <- DESeqDataSetFromMatrix(countData = FC_trx2wk_T, colData = info_trx2wk_T, design = ~Condition)
ddstrx <- DESeqDataSetFromMatrix(countData = FC_trx, colData = info_trx, design = ~Condition)
ddstrx2wk_all <- DESeqDataSetFromMatrix(countData = FC_trx2wk_all, colData = info_trx2wk_all, design = ~Condition)
ddstrx_all <- DESeqDataSetFromMatrix(countData = FC_trx_all, colData = info_trx_all, design = ~Condition)
ddsfig4 <- DESeqDataSetFromMatrix(countData = FC_fig_4, colData = info_fig4, design = ~Condition)


#remove lowly expressed genes
keepdoxo_all <- rowSums2(counts(ddsdoxo_all)) >= 150
ddsdoxo_all <- ddsdoxo_all[keepdoxo_all,]
keepdoxo_KO <- rowSums2(counts(ddsdoxo_KO)) >= 90
ddsdoxo_KO <- ddsdoxo_KO[keepdoxo_KO,]
keepdoxo_NT <- rowSums2(counts(ddsdoxo_NT)) >= 90
ddsdoxo_NT <- ddsdoxo_NT[keepdoxo_NT,]
keeptrx2wk_LSK <- rowSums2(counts(ddstrx2wk_LSK)) >= 120
ddstrx2wk_LSK <- ddstrx2wk_LSK[keeptrx2wk_LSK,]
keeptrx2wk_B <- rowSums2(counts(ddstrx2wk_B)) >= 90
ddstrx2wk_B <- ddstrx2wk_B[keeptrx2wk_B,]
keeptrx2wk_T <- rowSums2(counts(ddstrx2wk_T)) >= 90
ddstrx2wk_T <- ddstrx2wk_T[keeptrx2wk_T,]
keeptrx <- rowSums2(counts(ddstrx)) >= 90
ddstrx <- ddstrx[keeptrx,]
keeptrx2wk_all <- rowSums2(counts(ddstrx2wk_all)) >= 400
ddstrx2wk_all <- ddstrx2wk_all[keeptrx2wk_all,]
keeptrx_all <- rowSums2(counts(ddstrx_all)) >= 200
ddstrx_all <- ddstrx_all[keeptrx_all,]
keepfig4 <- rowSums2(counts(ddsfig4)) >= 400
ddsfig4 <- ddsfig4[keepfig4,]

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
DEtrx2wk_LSK <- DESeq(ddstrx2wk_LSK)
DEtrx2wk_LSKshr <- lfcShrink(DEtrx2wk_LSK, coef=2, type='apeglm')
plotDispEsts(DEtrx2wk_LSK, ylim =c(1e-4,1e2))
DEtrx2wk_B <- DESeq(ddstrx2wk_B)
DEtrx2wk_Bshr <- lfcShrink(DEtrx2wk_B, coef=2, type='apeglm')
plotDispEsts(DEtrx2wk_B, ylim =c(1e-4,1e2))
DEtrx2wk_T <- DESeq(ddstrx2wk_T)
DEtrx2wk_Tshr <- lfcShrink(DEtrx2wk_T, coef=2, type='apeglm')
plotDispEsts(DEtrx2wk_T, ylim =c(1e-4,1e2))
DEtrx <- DESeq(ddstrx)
DEtrxshr <- lfcShrink(DEtrx, coef=2, type='apeglm')
plotDispEsts(DEtrx, ylim =c(1e-4,1e2))
DEtrx2wk_all <- DESeq(ddstrx2wk_all)
DEtrx2wk_allshr <- lfcShrink(DEtrx2wk_all, coef=2, type='apeglm')
DEtrx_all <- DESeq(ddstrx_all)
DEtrx_allshr <- lfcShrink(DEtrx_all, coef=2, type='apeglm')
DEfig4 <- DESeq(ddsfig4)
DEfig4shr <- lfcShrink(DEfig4, coef=2, type='apeglm')

##MA plot for DEseq result
plotMA(results(DEdoxo_all, contrast = c('Condition','NT_het','NT_wt'), alpha = 0.05), ylim =c(-15,15),main="Doxoall: NT_het vs NT_wt")
plotMA(results(DEdoxo_all, contrast = c('Condition','Doxo_wt','NT_wt'), alpha = 0.05), ylim =c(-15,15),main="Doxoall: Doxo_wt vs NT_wt")
plotMA(results(DEdoxo_all, contrast = c('Condition','Doxo_het','NT_het'), alpha = 0.05), ylim =c(-15,15),main="Doxoall: Doxo_het vs NT_het")
plotMA(results(DEdoxo_all, contrast = c('Condition','Doxo_het','Doxo_wt'), alpha = 0.05), ylim =c(-15,15),main="Doxoall: Doxo_het vs Doxo_wt")
plotMA(results(DEdoxo_KO , contrast = c('Condition','Mut','WT'), alpha = 0.05), ylim =c(-15,15),main="DoxoKO: Mut vs WT")
plotMA(results(DEdoxo_NT , contrast = c('Condition','Mut','WT'), alpha = 0.05), ylim =c(-15,15),main="DoxoNT: Mut vs WT")
plotMA(results(DEtrx2wk_LSK, contrast = c('Condition','Mut','WT'), alpha = 0.05), ylim =c(-15,15),main="Trx2wk_LSK: Mut vs WT")
plotMA(results(DEtrx2wk_B, contrast = c('Condition','Mut','WT'), alpha = 0.05), ylim =c(-15,15),main="Trx2wk_B: Mut vs WT")
plotMA(results(DEtrx2wk_T, contrast = c('Condition','Mut','WT'), alpha = 0.05), ylim =c(-15,15),main="Trx2wk_T: Mut vs WT")
plotMA(results(DEtrx, contrast = c('Genotype','HET','WT'), alpha = 0.05), ylim =c(-15,15))
plotMA(results(DEtrx2wk_all, contrast = c('Condition','Mut_LSK','WT_LSK'), alpha = 0.05), ylim =c(-15,15),main="Trx2wk_LSK: Mut vs WT")
plotMA(results(DEtrx2wk_all, contrast = c('Condition','Mut_B','WT_B'), alpha = 0.05), ylim =c(-15,15),main="Trx2wk_LSK: Mut vs WT")
plotMA(results(DEtrx2wk_all, contrast = c('Condition','Mut_T','WT_T'), alpha = 0.05), ylim =c(-15,15),main="Trx2wk_LSK: Mut vs WT")


#export normalized read counts
normCounts_doxo_all <- counts(DEdoxo_all,normalized = T)
write.csv(normCounts_doxo_all,"~/Desktop/220624_rnaseq_Doxoall_cnt.csv")
Cnt_doxo_all <- read.csv("~/Desktop/220624_rnaseq_Doxoall_cnt.csv", row.names=1)
normCounts_doxo_KO <- counts(DEdoxo_KO,normalized = T)
write.csv(normCounts_doxo_KO,"~/Desktop/220624_rnaseq_DoxoKO_cnt.csv")
Cnt_doxo_KO <- read.csv("~/Desktop/220624_rnaseq_DoxoKO_cnt.csv", row.names=1)
normCounts_doxo_NT <- counts(DEdoxo_NT,normalized = T)
write.csv(normCounts_doxo_NT,"~/Desktop/220624_rnaseq_DoxoNT_cnt.csv")
Cnt_doxo_NT <- read.csv("~/Desktop/220624_rnaseq_DoxoNT_cnt.csv", row.names=1)
normCounts_trx2wk_LSK <- counts(DEtrx2wk_LSK,normalized = T)
write.csv(normCounts_trx2wk_LSK,"~/Desktop/220624_rnaseq_trx2wk_LSK_cnt.csv")
Cnt_trx2wk_LSK <- read.csv("~/Desktop/220624_rnaseq_trx2wk_LSK_cnt.csv", row.names=1)
normCounts_trx2wk_B <- counts(DEtrx2wk_B,normalized = T)
write.csv(normCounts_trx2wk_B,"~/Desktop/220624_rnaseq_trx2wk_B_cnt.csv")
Cnt_trx2wk_B <- read.csv("~/Desktop/220624_rnaseq_trx2wk_B_cnt.csv", row.names=1)
normCounts_trx2wk_T <- counts(DEtrx2wk_T,normalized = T)
write.csv(normCounts_trx2wk_T,"~/Desktop/220624_rnaseq_trx2wk_T_cnt.csv")
Cnt_trx2wk_T <- read.csv("~/Desktop/220624_rnaseq_trx2wk_T_cnt.csv", row.names=1)
normCounts_trx <- counts(DEtrx,normalized = T)
write.csv(normCounts_trx,"~/Desktop/220624_rnaseq_trx_cnt.csv")
Cnt_trx <- read.csv("~/Desktop/220624_rnaseq_trx_cnt.csv", row.names=1)
normCounts_trx2wk_all <- counts(DEtrx2wk_all,normalized = T)
write.csv(normCounts_trx2wk_all,"~/Desktop/220624_rnaseq_trx2wk_all_cnt.csv")
Cnt_trx2wk_all <- read.csv("~/Desktop/220624_rnaseq_trx2wk_all_cnt.csv", row.names=1)
normCounts_trx_all <- counts(DEtrx_all,normalized = T)
write.csv(normCounts_trx_all,"~/Desktop/220624_rnaseq_trx_all_cnt.csv")
Cnt_trx_all <- read.csv("~/Desktop/220624_rnaseq_trx_all_cnt.csv", row.names=1)
normCounts_fig4 <- counts(DEfig4,normalized = T)
write.csv(normCounts_fig4,"~/Desktop/220624_rnaseq_fig4_cnt.csv")
Cnt_fig4 <- read.csv("~/Desktop/220624_rnaseq_fig4_cnt.csv", row.names=1)

################
#PCA analysis
tnormCounts_doxo_all <- t(normCounts_doxo_all)
pca_doxo_all <- prcomp(tnormCounts_doxo_all, scale. = TRUE)
tnormCounts_doxo_KO <- t(normCounts_doxo_KO)
pca_doxo_KO <- prcomp(tnormCounts_doxo_KO, scale. = TRUE)
tnormCounts_doxo_NT <- t(normCounts_doxo_NT)
pca_doxo_NT <- prcomp(tnormCounts_doxo_NT, scale. = TRUE)
tnormCounts_trx2wk_LSK <- t(normCounts_trx2wk_LSK)
pca_trx2wk_LSK <- prcomp(tnormCounts_trx2wk_LSK, scale. = TRUE)
tnormCounts_trx2wk_B <- t(normCounts_trx2wk_B)
pca_trx2wk_B <- prcomp(tnormCounts_trx2wk_B, scale. = TRUE)
tnormCounts_trx2wk_T <- t(normCounts_trx2wk_T)
pca_trx2wk_T <- prcomp(tnormCounts_trx2wk_T, scale. = TRUE)
tnormCounts_trx <- t(normCounts_trx)
pca_trx <- prcomp(tnormCounts_trx, scale. = TRUE)
tnormCounts_trx2wk_all <- t(normCounts_trx2wk_all)
pca_trx2wk_all <- prcomp(tnormCounts_trx2wk_all, scale. = TRUE)
tnormCounts_trx_all <- t(normCounts_trx_all)
pca_trx_all <- prcomp(tnormCounts_trx_all, scale. = TRUE)
tnormCounts_fig4 <- t(normCounts_fig4)
pca_fig4 <- prcomp(tnormCounts_fig4, scale. = TRUE)
screeplot(pca_doxo_all, type = "l", main = "Screeplot for O2_KO")
abline(h=2000, col = 'red', lty =2)

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

summary(pca_doxo_KO)
pca_doxo_KO$x
pca_doxo_KO <- data.frame(info_doxo_KO[,2], pca_doxo_KO$x[,1:2])
colnames(pca_doxo_KO) <- c('Condition', 'PC1', 'PC2')
ggplot(pca_doxo_KO, aes(PC1, PC2, group=Condition)) + geom_point(size=8, aes(shape = Condition, color = Condition))+ 
  scale_shape_manual(values = c(15, 15))+
  scale_color_manual(values = c("red3","hotpink2"))+
  labs(x= "PC1", y= "PC2") + 
  theme(axis.text.y= element_text(size =16), axis.text.x= element_text(size =16), axis.title.x.bottom = element_text(size=20),axis.title.y = element_text(size=20))+
  xlim(-150,150) + ylim(-100, 100)

summary(pca_doxo_NT)
pca_doxo_NT$x
pca_doxo_NT <- data.frame(info_doxo_NT[,2], pca_doxo_NT$x[,1:2])
colnames(pca_doxo_NT) <- c('Condition', 'PC1', 'PC2')
ggplot(pca_doxo_NT, aes(PC1, PC2, group=Condition)) + geom_point(size=8, aes(shape = Condition, color = Condition))+ 
  scale_shape_manual(values = c(15, 15))+
  scale_color_manual(values = c("blue4","steelblue3"))+
  labs(x= "PC1", y= "PC2") + 
  theme(axis.text.y= element_text(size =16), axis.text.x= element_text(size =16), axis.title.x.bottom = element_text(size=20),axis.title.y = element_text(size=20))+
  xlim(-150,150) + ylim(-120, 120)

summary(pca_trx2wk_LSK)
pca_trx2wk_LSK$x
pca_trx2wk_LSK <- data.frame(info_trx2wk_LSK[,2], pca_trx2wk_LSK$x[,1:2])
colnames(pca_trx2wk_LSK) <- c('Condition', 'PC1', 'PC2')
ggplot(pca_trx2wk_LSK, aes(PC1, PC2, group=Condition)) + geom_point(size=8, aes(shape = Condition, color = Condition))+ 
  scale_shape_manual(values = c(15,16,15,16))+
  scale_color_manual(values = c("blue4","steelblue3","hotpink2","red3"))+
  labs(x= "PC1", y= "PC2") + 
  theme(axis.text.y= element_text(size =16), axis.text.x= element_text(size =16), axis.title.x.bottom = element_text(size=20),axis.title.y = element_text(size=20))+
  xlim(-150,150) + ylim(-150, 150)

summary(pca_trx2wk_B)
pca_trx2wk_B$x
pca_trx2wk_B <- data.frame(info_trx2wk_B[,2], pca_trx2wk_B$x[,1:2])
colnames(pca_trx2wk_B) <- c('Condition', 'PC1', 'PC2')
ggplot(pca_trx2wk_B, aes(PC1, PC2, group=Condition)) + geom_point(size=8, aes(shape = Condition, color = Condition))+ 
  scale_shape_manual(values = c(15, 16))+
  scale_color_manual(values = c("blue4","red3"))+
  labs(x= "PC1", y= "PC2") + 
  theme(axis.text.y= element_text(size =16), axis.text.x= element_text(size =16), axis.title.x.bottom = element_text(size=20),axis.title.y = element_text(size=20))+
  xlim(-150,150) + ylim(-100, 100)

summary(pca_trx2wk_T)
pca_trx2wk_T$x
pca_trx2wk_T <- data.frame(info_trx2wk_T[,2], pca_trx2wk_T$x[,1:2])
colnames(pca_trx2wk_T) <- c('Condition', 'PC1', 'PC2')
ggplot(pca_trx2wk_T, aes(PC1, PC2, group=Condition)) + geom_point(size=8, aes(shape = Condition, color = Condition))+ 
  scale_shape_manual(values = c(15, 16))+
  scale_color_manual(values = c("blue4","red3"))+
  labs(x= "PC1", y= "PC2") + 
  theme(axis.text.y= element_text(size =16), axis.text.x= element_text(size =16), axis.title.x.bottom = element_text(size=20),axis.title.y = element_text(size=20))+
  xlim(-150,150) + ylim(-100, 100)

summary(pca_trx2wk_all)
pca_trx2wk_all$x
pca_trx2wk_all <- data.frame(info_trx2wk_all[,2], pca_trx2wk_all$x[,1:2])
colnames(pca_trx2wk_all) <- c('Condition', 'PC1', 'PC2')
ggplot(pca_trx2wk_all, aes(PC1, PC2, group=Condition)) + geom_point(size=4, aes(shape = Condition, color = Condition))+ 
  scale_shape_manual(values = c(15, 15,15,16,16,16))+
  scale_color_manual(values = c("red3","gray25","blue4","hotpink2","gray47","steelblue3"))+
  labs(x= "PC1", y= "PC2") + 
  theme(axis.text.y= element_text(size =16), axis.text.x= element_text(size =16), axis.title.x.bottom = element_text(size=20),axis.title.y = element_text(size=20))+
  xlim(-110,110) + ylim(-120, 80)

summary(pca_trx)
pca_trx$x
pca_trx1 <- data.frame(infotrx1[,2], pca_trx$x[,1:2])
colnames(pca_trx1) <- c('Treatment', 'PC1', 'PC2')
ggplot(pca_trx1, aes(PC1, PC2, color = Treatment)) + geom_point(size=8) + 
  labs(x= "PC1", y= "PC2") + 
  theme(axis.text.y= element_text(size =16), axis.text.x= element_text(size =16), 
        axis.title.x.bottom = element_text(size=20),axis.title.y = element_text(size=20))+
  xlim(-150,150) + ylim(-150, 150)

summary(pca_trx_all)
pca_trx_all$x
pca_trx_all <- data.frame(info_trx_all[,2], pca_trx_all$x[,1:2])
colnames(pca_trx_all) <- c('Condition', 'PC1', 'PC2')
ggplot(pca_trx_all, aes(PC1, PC2, group=Condition)) + geom_point(size=8, aes(shape = Condition, color = Condition))+ 
  scale_shape_manual(values = c(15, 16,15,16))+
  scale_color_manual(values = c("red3","hotpink2", "blue4","steelblue3","blue4"))+
  labs(x= "PC1", y= "PC2") + 
  theme(axis.text.y= element_text(size =16), axis.text.x= element_text(size =16), axis.title.x.bottom = element_text(size=20),axis.title.y = element_text(size=20))+
  xlim(-150,150) + ylim(-120, 120)
  
pca_fig4$x
pca_fig4 <- data.frame(info_fig4[,2], pca_fig4$x[,1:2])
colnames(pca_fig4) <- c('Condition', 'PC1', 'PC2')
ggplot(pca_fig4, aes(PC1, PC2, group=Condition)) + geom_point(size=5, aes(shape = Condition, color = Condition))+ 
  scale_shape_manual(values = c(15, 15,15,16,16,16))+
  scale_color_manual(values = c("gray34","gray49","hotpink2","red3", "blue4","steelblue3"))+
  labs(x= "PC1", y= "PC2") + 
  theme(axis.text.y= element_text(size =16), axis.text.x= element_text(size =16), axis.title.x.bottom = element_text(size=20),axis.title.y = element_text(size=20))+
  xlim(-100,150) + ylim(-80, 80)

################
#DESeq results
res_doxo_all_doxo <- results(DEdoxo_all, contrast = c("Condition", "Doxo_het", "Doxo_wt"), alpha = 0.05, cooksCutoff = FALSE)
res_doxo_all_NT <- results(DEdoxo_all, contrast = c("Condition", "NT_het", "NT_wt"), alpha = 0.05, cooksCutoff = FALSE)
res_doxo_KO <- results(DEdoxo_KO, contrast = c("Condition", "Mut", "WT"), alpha = 0.05, cooksCutoff = FALSE)
res_doxo_NT <- results(DEdoxo_NT, contrast = c("Condition", "Mut", "WT"), alpha = 0.05, cooksCutoff = FALSE)
res_trx2wk_LSK <- results(DEtrx2wk_LSK , contrast = c("Condition", "Mut", "WT"), alpha = 0.05, cooksCutoff = FALSE)
res_trx2wk_B <- results(DEtrx2wk_B , contrast = c("Condition", "Mut", "WT"), alpha = 0.05, cooksCutoff = FALSE)
res_trx2wk_T <- results(DEtrx2wk_T , contrast = c("Condition", "Mut", "WT"), alpha = 0.05, cooksCutoff = FALSE)
res_trx <- results(DEtrx, contrast = c("Condition","Mut", "WT"), alpha = 0.05, cooksCutoff = FALSE)
res_trx2wk_all_LSK <- results(DEtrx2wk_all , contrast = c("Condition", "Mut_LSK", "WT_LSK"), alpha = 0.05, cooksCutoff = FALSE)
res_trx2wk_all_B <- results(DEtrx2wk_all , contrast = c("Condition", "Mut_B", "WT_B"), alpha = 0.05, cooksCutoff = FALSE)
res_trx2wk_all_T <- results(DEtrx2wk_all , contrast = c("Condition", "Mut_T", "WT_T"), alpha = 0.05, cooksCutoff = FALSE)
res_trx_all_2wk <- results(DEtrx_all, contrast = c("Condition","Mut_2wk","WT_2wk" ), alpha = 0.05, cooksCutoff = FALSE)
res_trx_all_5mo <- results(DEtrx_all, contrast = c("Condition","Mut_5mo","WT_5mo" ), alpha = 0.05, cooksCutoff = FALSE)

#output DESeq results
res_doxo_all_doxo <- res_doxo_all_doxo[order(res_doxo_all_doxo$padj),]
res_doxo_all_NT <- res_doxo_all_NT[order(res_doxo_all_NT$padj),]
res_doxo_KO <- res_doxo_KO[order(res_doxo_KO$padj),]
res_doxo_NT <- res_doxo_NT[order(res_doxo_NT$padj),]
res_trx2wk_LSK <- res_trx2wk_LSK[order(res_trx2wk_LSK$padj),]
res_trx2wk_B <- res_trx2wk_B[order(res_trx2wk_B$padj),]
res_trx2wk_T <- res_trx2wk_T[order(res_trx2wk_T$padj),]
res_trx <- res_trx[order(res_trx$padj),]
res_trx2wk_all_LSK <- res_trx2wk_all_LSK[order(res_trx2wk_all_LSK$padj),]
res_trx2wk_all_B <- res_trx2wk_all_B[order(res_trx2wk_all_B$padj),]
res_trx2wk_all_T <- res_trx2wk_all_T[order(res_trx2wk_all_T$padj),]
res_trx_all_2wk <- res_trx_all_2wk[order(res_trx_all_2wk$padj),]
res_trx_all_5mo <- res_trx_all_5mo[order(res_trx_all_5mo$padj),]

write.csv(res_doxo_all_doxo, "~/Desktop/220624_rnaseq_doxoall_doxo_pval.csv")
write.csv(res_doxo_all_NT, "~/Desktop/220624_rnaseq_doxoall_NT_pval.csv")
write.csv(res_doxo_KO, "~/Desktop/220624_rnaseq_doxo_KO_pval.csv")
write.csv(res_doxo_NT, "~/Desktop/220624_rnaseq_doxo_NT_pval.csv")
write.csv(res_trx2wk_LSK, "~/Desktop/220624_rnaseq_trx2wk_LSK_pval.csv")
write.csv(res_trx2wk_B, "~/Desktop/220624_rnaseq_trx2wk_B_pval.csv")
write.csv(res_trx2wk_T, "~/Desktop/220624_rnaseq_trx2wk_T_pval.csv")
write.csv(res_trx, "~/Desktop/220624_rnaseq_trx_pval.csv")
write.csv(res_trx2wk_all_LSK, "~/Desktop/220624_rnaseq_trx2wk_all_LSK_pval.csv")
write.csv(res_trx2wk_all_B, "~/Desktop/220624_rnaseq_trx2wk_all_B_pval.csv")
write.csv(res_trx2wk_all_T, "~/Desktop/220624_rnaseq_trx2wk_all_T_pval.csv")
write.csv(res_trx_all_2wk, "~/Desktop/220624_rnaseq_trx_all_2wk_pval.csv")
write.csv(res_trx_all_5mo, "~/Desktop/220624_rnaseq_trx_all_5mo_pval.csv")


Pval_all_doxo <- read.csv("~/Desktop/220624_rnaseq_doxoall_doxo_pval.csv", row.names = 1)
Pval_all_NT <- read.csv("~/Desktop/220624_rnaseq_doxoall_NT_pval.csv", row.names = 1)
Pval_doxo_KO <- read.csv("~/Desktop/220624_rnaseq_doxo_KO_pval.csv", row.names = 1)
Pval_doxo_NT <- read.csv("~/Desktop/220624_rnaseq_doxo_NT_pval.csv", row.names = 1)
Pval_trx2wk_LSK <- read.csv("~/Desktop/220624_rnaseq_trx2wk_LSK_pval.csv", row.names = 1)
Pval_trx2wk_B <- read.csv("~/Desktop/220624_rnaseq_trx2wk_B_pval.csv", row.names = 1)
Pval_trx2wk_T <- read.csv("~/Desktop/220624_rnaseq_trx2wk_T_pval.csv", row.names = 1)
Pval_trx <- read.csv("~/Desktop/220624_rnaseq_trx_pval.csv", row.names=1)
Pval_trx2wk_all_LSK <- read.csv("~/Desktop/220624_rnaseq_trx2wk_all_LSK_pval.csv", row.names = 1)
Pval_trx2wk_all_B <- read.csv("~/Desktop/220624_rnaseq_trx2wk_all_B_pval.csv", row.names = 1)
Pval_trx2wk_all_T <- read.csv("~/Desktop/220624_rnaseq_trx2wk_all_T_pval.csv", row.names = 1)
Pval_trx_all_2wk <- read.csv("~/Desktop/220624_rnaseq_trx_all_2wk_pval.csv", row.names=1)
Pval_trx_all_5mo <- read.csv("~/Desktop/220624_rnaseq_trx_all_5mo_pval.csv", row.names=1)

##Reoragnize DESeq result
Pval_all_doxo$sig <- ifelse(Pval_all_doxo$padj <= 0.05, "yes", "no")
Pval_all_NT$sig <- ifelse(Pval_all_NT$padj <= 0.05, "yes", "no")
Pval_doxo_KO$sig <- ifelse(Pval_doxo_KO$padj <= 0.05, "yes", "no")
Pval_doxo_NT$sig <- ifelse(Pval_doxo_NT$padj <= 0.05, "yes", "no")
Pval_trx2wk_LSK$sig <- ifelse(Pval_trx2wk_LSK$padj <= 0.05, "yes", "no")
Pval_trx2wk_B$sig <- ifelse(Pval_trx2wk_B$padj <= 0.05, "yes", "no")
Pval_trx2wk_T$sig <- ifelse(Pval_trx2wk_T$padj <= 0.05, "yes", "no")
Pval_trx$sig <- ifelse(Pval_trx$padj <= 0.05, "yes", "no")
Pval_trx2wk_all_LSK$sig <- ifelse(Pval_trx2wk_all_LSK$padj <= 0.05, "yes", "no")
Pval_trx2wk_all_B$sig <- ifelse(Pval_trx2wk_all_B$padj <= 0.05, "yes", "no")
Pval_trx2wk_all_T$sig <- ifelse(Pval_trx2wk_all_T$padj <= 0.05, "yes", "no")
Pval_trx_all_2wk$sig <- ifelse(Pval_trx_all_2wk$padj <= 0.05, "yes", "no")
Pval_trx_all_5mo$sig <- ifelse(Pval_trx_all_5mo$padj <= 0.05, "yes", "no")
Pval_S$sig <- ifelse(Pval_S$padj <= 0.05, "yes", "no")


##Filter outliers
Pval_all_doxo<- na.omit(Pval_all_doxo)
Pval_all_NT<- na.omit(Pval_all_NT)
Pval_doxo_KO<- na.omit(Pval_doxo_KO)
Pval_doxo_NT<- na.omit(Pval_doxo_NT)
Pval_trx2wk_LSK<- na.omit(Pval_trx2wk_LSK)
Pval_trx2wk_B<- na.omit(Pval_trx2wk_B)
Pval_trx2wk_T<- na.omit(Pval_trx2wk_T)
Pval_trx<- na.omit(Pval_trx)
Pval_trx2wk_all_LSK<- na.omit(Pval_trx2wk_all_LSK)
Pval_trx2wk_all_B<- na.omit(Pval_trx2wk_all_B)
Pval_trx2wk_all_T<- na.omit(Pval_trx2wk_all_T)
Pval_trx_all_2wk<- na.omit(Pval_trx_all_2wk)
Pval_trx_all_5mo<- na.omit(Pval_trx_all_5mo)
Pval_S<- na.omit(Pval_S)

#remove fold change outlier
Pval_all_doxo <- filter(Pval_all_doxo,log2FoldChange <10, -log2FoldChange<10)
Pval_all_NT <- filter(Pval_all_NT,log2FoldChange <10, -log2FoldChange<10)
Pval_doxo_KO <- filter(Pval_doxo_KO,log2FoldChange <10, -log2FoldChange<10)
Pval_doxo_NT <- filter(Pval_doxo_NT,log2FoldChange <10, -log2FoldChange<10)
Pval_trx2wk_LSK <- filter(Pval_trx2wk_LSK,log2FoldChange <8, -log2FoldChange<8)
Pval_trx2wk_B <- filter(Pval_trx2wk_B,log2FoldChange <8, -log2FoldChange<8)
Pval_trx2wk_T <- filter(Pval_trx2wk_T,log2FoldChange <8, -log2FoldChange<8)
Pval_trx <- filter(Pval_trx, log2FoldChange<7, -log2FoldChange<7)
Pval_trx2wk_all_LSK <- filter(Pval_trx2wk_all_LSK, log2FoldChange <8, -log2FoldChange<8)
Pval_trx2wk_all_B <- filter(Pval_trx2wk_all_B,log2FoldChange <8, -log2FoldChange<8)
Pval_trx2wk_all_T <- filter(Pval_trx2wk_all_T,log2FoldChange <8, -log2FoldChange<8)
Pval_trx_all_2wk <- filter(Pval_trx_all_2wk, log2FoldChange<10, -log2FoldChange<10)
Pval_trx_all_5mo <- filter(Pval_trx_all_5mo, log2FoldChange<10, -log2FoldChange<10)
Pval_S <- filter(Pval_S, log2FoldChange <7, -log2FoldChange<7)

##Esembl id to gene name
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
Pval_trx2wk_LSK <- merge(Pval_trx2wk_LSK, Ensembl.id.list ,by=0)
row.names(Pval_trx2wk_LSK) <- Pval_trx2wk_LSK[,9]
Pval_trx2wk_LSK <- Pval_trx2wk_LSK[,-9:-10]
Pval_trx2wk_LSK <- Pval_trx2wk_LSK[,-1]
Pval_trx2wk_B <- merge(Pval_trx2wk_B, Ensembl.id.list ,by=0)
row.names(Pval_trx2wk_B) <- Pval_trx2wk_B[,9]
Pval_trx2wk_B <- Pval_trx2wk_B[,-9:-10]
Pval_trx2wk_B <- Pval_trx2wk_B[,-1]
Pval_trx2wk_T <- merge(Pval_trx2wk_T, Ensembl.id.list ,by=0)
row.names(Pval_trx2wk_T) <- Pval_trx2wk_T[,9]
Pval_trx2wk_T <- Pval_trx2wk_T[,-9:-10]
Pval_trx2wk_T <- Pval_trx2wk_T[,-1]
Pval_trx <- merge(Pval_trx, Ensembl.id.list ,by=0)
row.names(Pval_trx) <- Pval_trx[,9]
Pval_trx <- Pval_trx[,-9:-10]

Pval_trx2wk_all_LSK <- merge(Pval_trx2wk_all_LSK, Ensembl.id.list ,by=0)
row.names(Pval_trx2wk_all_LSK) <- Pval_trx2wk_all_LSK[,9]
Pval_trx2wk_all_LSK <- Pval_trx2wk_all_LSK[,-9:-10]
Pval_trx2wk_all_LSK <- Pval_trx2wk_all_LSK[,-1]
Pval_trx2wk_all_B <- merge(Pval_trx2wk_all_B, Ensembl.id.list ,by=0)
row.names(Pval_trx2wk_all_B) <- Pval_trx2wk_all_B[,9]
Pval_trx2wk_all_B <- Pval_trx2wk_all_B[,-9:-10]
Pval_trx2wk_all_B <- Pval_trx2wk_all_B[,-1]
Pval_trx2wk_all_T <- merge(Pval_trx2wk_all_T, Ensembl.id.list ,by=0)
row.names(Pval_trx2wk_all_T) <- Pval_trx2wk_all_T[,9]
Pval_trx2wk_all_T <- Pval_trx2wk_all_T[,-9:-10]
Pval_trx2wk_all_T <- Pval_trx2wk_all_T[,-1]
                 
Pval_trx_all_2wk <- merge(Pval_trx_all_2wk, Ensembl.id.list ,by=0)
row.names(Pval_trx_all_2wk) <- Pval_trx_all_2wk[,9]
Pval_trx_all_2wk <- Pval_trx_all_2wk[,-9:-10]
Pval_trx_all_5mo <- merge(Pval_trx_all_5mo, Ensembl.id.list ,by=0)
row.names(Pval_trx_all_5mo) <- Pval_trx_all_5mo[,9]
Pval_trx_all_5mo <- Pval_trx_all_5mo[,-9:-10]

write.csv(res_doxo_all_doxo, "~/Desktop/220624_rnaseq_doxoall_doxo_pval.csv")
write.csv(res_doxo_all_NT, "~/Desktop/220624_rnaseq_doxoall_NT_pval.csv")

write.csv(Pval_all_doxo, "~/Desktop/Pval_all_doxo.csv")
write.csv(Pval_all_NT, "~/Desktop/Pval_all_NT.csv")
write.csv(Pval_doxo_KO, "~/Desktop/Pval_doxo_KO.csv")
write.csv(Pval_doxo_NT, "~/Desktop/Pval_doxo_NT.csv")
write.csv(Pval_trx2wk_LSK, "~/Desktop/Pval_trx2wk_LSK.csv")
write.csv(Pval_trx2wk_B, "~/Desktop/Pval_trx2wk_B.csv")
write.csv(Pval_trx2wk_T, "~/Desktop/Pval_trx2wk_T.csv")
write.csv(Pval_trx, "~/Desktop/Pval_trx.csv")
write.csv(Pval_trx2wk_all_LSK, "~/Desktop/Pval_trx2wk_all_LSK.csv")
write.csv(Pval_trx2wk_all_B, "~/Desktop/Pval_trx2wk_all_B.csv")
write.csv(Pval_trx2wk_all_T, "~/Desktop/Pval_trx2wk_all_T.csv")
write.csv(Pval_trx_all_2wk, "~/Desktop/Pval_trx_all_2wk.csv")
write.csv(Pval_trx_all_5mo, "~/Desktop/Pval_trx_all_5mo.csv")

Pval_all_doxo <- read.csv("~/Desktop/Pval_all_doxo.csv", row.names=1)
Pval_all_NT <- read.csv("~/Desktop/Pval_all_NT.csv", row.names=1)
Pval_doxo_KO <- read.csv("~/Desktop/Pval_doxo_KO.csv", row.names=1)
Pval_doxo_NT <- read.csv("~/Desktop/Pval_doxo_NT.csv", row.names=1)
Pval_trx2wk_LSK <- read.csv("~/Desktop/Pval_trx2wk_LSK.csv", row.names=1)
Pval_trx2wk_B <- read.csv("~/Desktop/Pval_trx2wk_B.csv", row.names=1)
Pval_trx2wk_T <- read.csv("~/Desktop/Pval_trx2wk_T.csv", row.names=1)
Pval_trx <- read.csv("~/Desktop/Pval_trx.csv", row.names=1)
Pval_trx2wk_all_LSK <- read.csv("~/Desktop/Pval_trx2wk_all_LSK.csv", row.names=1)
Pval_trx2wk_all_B <- read.csv("~/Desktop/Pval_trx2wk_all_B.csv", row.names=1)
Pval_trx2wk_all_T <- read.csv("~/Desktop/Pval_trx2wk_all_T.csv", row.names=1)
Pval_trx_all_2wk <- read.csv("~/Desktop/Pval_trx_all_2wk.csv", row.names=1)
Pval_trx_all_5mo <- read.csv("~/Desktop/Pval_trx_all_5mo.csv", row.names=1)

Pval_all_doxo <- Pval_all_doxo[order(Pval_all_doxo$padj),]
Pval_all_NT <- Pval_all_NT[order(Pval_all_NT$padj),]
Pval_doxo_KO <- Pval_doxo_KO[order(Pval_doxo_KO$padj),]
Pval_doxo_NT <- Pval_doxo_NT[order(Pval_doxo_NT$padj),]
Pval_trx2wk_LSK <- Pval_trx2wk_LSK[order(Pval_trx2wk_LSK$padj),]
Pval_trx2wk_B <- Pval_trx2wk_B[order(Pval_trx2wk_B$padj),]
Pval_trx2wk_T <- Pval_trx2wk_T[order(Pval_trx2wk_T$padj),]
Pval_trx <- Pval_trx[order(Pval_trx$padj),]
Pval_trx2wk_all_LSK <- Pval_trx2wk_all_LSK[order(Pval_trx2wk_all_LSK$padj),]
Pval_trx2wk_all_B <- Pval_trx2wk_all_B[order(Pval_trx2wk_all_B$padj),]
Pval_trx2wk_all_T <- Pval_trx2wk_all_T[order(Pval_trx2wk_all_T$padj),]
Pval_trx_all_2wk <- Pval_trx_all_2wk[order(Pval_trx_all_2wk$padj),]
Pval_trx_all_5mo <- Pval_trx_all_5mo[order(Pval_trx_all_5mo$padj),]

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
EnhancedVolcano(Pval_doxo_KO, lab = NA, x='log2FoldChange',y='padj',
                title = 'Doxorubicin: Mut vs WT', pCutoff = 5e-2, pointSize = 1.5,
                col=c('darkgrey','darkgrey', 'darkgrey','red'), colAlpha = 0.5,
                FCcutoff = 1.25, legendPosition = 'right') + xlim(-10,10) + ylim(0,15)
EnhancedVolcano(Pval_doxo_NT, lab = NA, x='log2FoldChange',y='padj',
                title = 'Germline: Mut vs WT', pCutoff = 5e-2, pointSize = 1.5,
                col=c('darkgrey','darkgrey', 'darkgrey','red'), colAlpha = 0.5,
                FCcutoff = 1.25, legendPosition = 'right') + xlim(-10,10) + ylim(0,10)
EnhancedVolcano(Pval_trx2wk_LSK, lab = NA, x='log2FoldChange',y='padj',
                title = 'Trx2wk_LSK: Mut vs WT', pCutoff = 5e-2, pointSize = 1.5,
                col=c('darkgrey','darkgrey', 'darkgrey','red'), colAlpha = 0.5,
                FCcutoff = 1.25, legendPosition = 'right') + xlim(-10,10) + ylim(0,10)
EnhancedVolcano(Pval_trx2wk_B, lab = NA, x='log2FoldChange',y='padj',
                title = 'Trx2wk_Bcell: Mut vs WT', pCutoff = 5e-2, pointSize = 1.5,
                col=c('darkgrey','darkgrey', 'darkgrey','red'), colAlpha = 0.5,
                FCcutoff = 1.25, legendPosition = 'right') + xlim(-10,10) + ylim(0,10)
EnhancedVolcano(Pval_trx2wk_T, lab = NA, x='log2FoldChange',y='padj',
                title = 'Trx2wk_Tcell: Mut vs WT', pCutoff = 5e-2, pointSize = 1.5,
                col=c('darkgrey','darkgrey', 'darkgrey','red'), colAlpha = 0.5,
                FCcutoff = 1.25, legendPosition = 'right') + xlim(-10,10) + ylim(0,10)
EnhancedVolcano(Pval_trx_all_2wk, lab = NA, x='log2FoldChange',y='padj',
                title = 'Trx_2wk: Mut vs WT', pCutoff = 5e-2, pointSize = 1.5,
                col=c('darkgrey','darkgrey', 'darkgrey','red'), colAlpha = 0.5,
                FCcutoff = 1.25, legendPosition = 'right') + xlim(-10,10) + ylim(0,10)
EnhancedVolcano(Pval_trx_all_5mo, lab = NA, x='log2FoldChange',y='padj',
                title = 'Trx_5mo: Mut vs WT', pCutoff = 5e-2, pointSize = 1.5,
                col=c('darkgrey','darkgrey', 'darkgrey','red'), colAlpha = 0.5,
                FCcutoff = 1.25, legendPosition = 'right') + xlim(-10,10) + ylim(0,10)
EnhancedVolcano(Pval_trx, lab = NA, x='log2FoldChange',y='padj',
                title = 'Trx: Mut vs WT', pCutoff = 5e-2, pointSize = 1.5,
                col=c('darkgrey','darkgrey', 'darkgrey','darkblue'), colAlpha = 0.5,
                FCcutoff = 1.25, legendPosition = 'right') + xlim(-10,10) + ylim(0,8)


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

Pval_trx2wk_LSK_f <- subset(Pval_trx2wk_LSK, baseMean >=150)
Genelist_trx2wk_LSK <- Pval_trx2wk_LSK_f$log2FoldChange
names(Genelist_trx2wk_LSK) <- row.names(Pval_trx2wk_LSK_f)
Genelist_trx2wk_LSK <- na.omit(Genelist_trx2wk_LSK)
Genelist_trx2wk_LSK <- sort(Genelist_trx2wk_LSK, decreasing = TRUE)
Gse_trx2wk_LSK <- gseGO(geneList = Genelist_trx2wk_LSK, ont = "BP", keyType = "ENSEMBL", 
                        minGSSize = 80, maxGSSize = 400, pvalueCutoff = 1, verbose = TRUE,
                        OrgDb = "org.Mm.eg.db",pAdjustMethod = "BH", eps=0)
Gsea_trx2wk_LSK <- as.data.frame(Gse_trx2wk_LSK)
write.csv(Gsea_trx2wk_LSK, "~/Desktop/Gsea_trx2wk_LSK.csv")
require(DOSE)
dotplot(Gse_trx2wk_LSK, font.size=8, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("GSEA: Transplantation_2wk")

Pval_trx2wk_B_f <- subset(Pval_trx2wk_B, baseMean >=150)
Genelist_trx2wk_B <- Pval_trx2wk_B_f$log2FoldChange
names(Genelist_trx2wk_B) <- row.names(Pval_trx2wk_B_f)
Genelist_trx2wk_B <- na.omit(Genelist_trx2wk_B)
Genelist_trx2wk_B <- sort(Genelist_trx2wk_B, decreasing = TRUE)
Gse_trx2wk_B <- gseGO(geneList = Genelist_trx2wk_B, ont = "ALL", keyType = "ENSEMBL", 
                      minGSSize = 20, maxGSSize = 300, pvalueCutoff = 0.05, verbose = TRUE,
                      OrgDb = "org.Mm.eg.db",pAdjustMethod = "BH", eps=0)
Gsea_trx2wk_B <- as.data.frame(Gse_trx2wk_B)
require(DOSE)
dotplot(Gse_trx2wk_B, font.size=8, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("GSEA: Doxorubicin Treatment")

Pval_trx_filter <- subset(Pval_trx, baseMean >=150)
Genelist_trx <- Pval_trx_filter$log2FoldChange
names(Genelist_trx) <- row.names(Pval_trx_filter)
Genelist_trx <- na.omit(Genelist_trx)
Genelist_trx <- sort(Genelist_trx, decreasing = TRUE)
Gse_trx <- gseGO(geneList = Genelist_trx, ont = "BP", keyType = "ENSEMBL", 
                      minGSSize = 80, maxGSSize = 400, pvalueCutoff = 0.05, verbose = TRUE,
                      OrgDb = "org.Mm.eg.db",pAdjustMethod = "BH", eps=0)
require(DOSE)
dotplot(Gse_trx, font.size=8, showCategory=10, split=".sign") + facet_grid(.~.sign) + ggtitle("GSEA: Transplantation")

Pval_trx2wk_all_LSK_f <- subset(Pval_trx2wk_all_LSK, baseMean >=150)
Genelist_trx2wk_all_LSK <- Pval_trx2wk_all_LSK_f$log2FoldChange
names(Genelist_trx2wk_all_LSK) <- row.names(Pval_trx2wk_all_LSK_f)
Genelist_trx2wk_all_LSK <- na.omit(Genelist_trx2wk_all_LSK)
Genelist_trx2wk_all_LSK <- sort(Genelist_trx2wk_all_LSK, decreasing = TRUE)
Gse_trx2wk_all_LSK <- gseGO(geneList = Genelist_trx2wk_all_LSK, ont = "BP", keyType = "ENSEMBL", 
                            minGSSize = 80, maxGSSize = 500, verbose = TRUE,
                            OrgDb = "org.Mm.eg.db", pvalueCutoff = 1 ,pAdjustMethod = "none", eps=0)
Gsea_trx2wk_all_LSK <- as.data.frame(Gse_trx2wk_all_LSK)

Pval_trx2wk_all_B_f <- subset(Pval_trx2wk_all_B, baseMean >=150)
Genelist_trx2wk_all_B <- Pval_trx2wk_all_B_f$log2FoldChange
names(Genelist_trx2wk_all_B) <- row.names(Pval_trx2wk_all_B_f)
Genelist_trx2wk_all_B <- na.omit(Genelist_trx2wk_all_B)
Genelist_trx2wk_all_B <- sort(Genelist_trx2wk_all_B, decreasing = TRUE)
Gse_trx2wk_all_B <- gseGO(geneList = Genelist_trx2wk_all_B, ont = "BP", keyType = "ENSEMBL", 
                            minGSSize = 80, maxGSSize = 400, pvalueCutoff = 0.05, verbose = TRUE,
                            OrgDb = "org.Mm.eg.db",pAdjustMethod = "BH", eps=0)
Gsea_trx2wk_all_B <- as.data.frame(Gse_trx2wk_all_B)

Pval_trx2wk_all_T_f <- subset(Pval_trx2wk_all_T, baseMean >=150)
Genelist_trx2wk_all_T <- Pval_trx2wk_all_T_f$log2FoldChange
names(Genelist_trx2wk_all_T) <- row.names(Pval_trx2wk_all_T_f)
Genelist_trx2wk_all_T <- na.omit(Genelist_trx2wk_all_T)
Genelist_trx2wk_all_T <- sort(Genelist_trx2wk_all_T, decreasing = TRUE)
Gse_trx2wk_all_T <- gseGO(geneList = Genelist_trx2wk_all_T, ont = "BP", keyType = "ENSEMBL", 
                          minGSSize = 80, maxGSSize = 400, pvalueCutoff = 0.05, verbose = TRUE,
                          OrgDb = "org.Mm.eg.db",pAdjustMethod = "BH", eps=0)
Gsea_trx2wk_all_T <- as.data.frame(Gse_trx2wk_all_T)

Pval_trx_all_2wk_f <- subset(Pval_trx_all_2wk, baseMean >=150)
Genelist_trx_all_2wk <- Pval_trx_all_2wk_f$log2FoldChange
names(Genelist_trx_all_2wk) <- row.names(Pval_trx_all_2wk_f)
Genelist_trx_all_2wk <- na.omit(Genelist_trx_all_2wk)
Genelist_trx_all_2wk <- sort(Genelist_trx_all_2wk, decreasing = TRUE)
Gse_trx_all_2wk <- gseGO(geneList = Genelist_trx_all_2wk, ont = "BP", keyType = "ENSEMBL", 
                         minGSSize = 80, maxGSSize = 500, verbose = TRUE,
                         OrgDb = "org.Mm.eg.db", pvalueCutoff = 1 ,pAdjustMethod = "none", eps=0)

Gsea_trx2wk_all_LSK <- as.data.frame(Gse_trx2wk_all_LSK)

Pval_trx_all_5mo_f <- subset(Pval_trx_all_5mo, baseMean >=150)
Genelist_trx_all_5mo <- Pval_trx_all_5mo_f$log2FoldChange
names(Genelist_trx_all_5mo) <- row.names(Pval_trx_all_5mo_f)
Genelist_trx_all_5mo <- na.omit(Genelist_trx_all_5mo)
Genelist_trx_all_5mo <- sort(Genelist_trx_all_5mo, decreasing = TRUE)
Gse_trx_all_5mo <- gseGO(geneList = Genelist_trx_all_5mo, ont = "BP", keyType = "ENSEMBL", 
                         minGSSize = 80, maxGSSize = 500, verbose = TRUE,
                         OrgDb = "org.Mm.eg.db", pvalueCutoff = 1 ,pAdjustMethod = "none", eps=0)
Gsea_trx2wk_all_LSK <- as.data.frame(Gse_trx2wk_all_LSK)

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
gseaplot2(Gse_trx2wk_LSK, geneSetID = 41, title = Gse_trx2wk_LSK$Description[41])
Gse_trx2wk_LSK$NES[41]
Gse_trx2wk_LSK$p.adjust[41]
gseaplot2(Gse_trx2wk_all_LSK , geneSetID = 244, title = Gse_trx2wk_all_LSK$Description[244])
Gse_trx2wk_all_LSK$NES[244]
Gse_trx2wk_all_LSK$p.adjust[244]
gseaplot2(Gse_trx, geneSetID = 49, title = Gse_trx$Description[49])
Gse_trx$NES[49]
gseaplot2(Gse_trx_all_2wk , geneSetID = 67, title = Gse_trx_all_2wk$Description[67])
Gse_trx_all_2wk$NES[67]
Gse_trx_all_2wk$p.adjust[67]
gseaplot2(Gse_trx_all_5mo , geneSetID = 369, title = Gse_trx_all_5mo$Description[369])
Gse_trx_all_5mo$NES[369]
Gse_trx_all_5mo$p.adjust[369]
gseaplot2(Gse_trx_all_5mo , geneSetID = 37, title = Gse_trx_all_5mo$Description[37])
Gse_trx_all_5mo$NES[37]
Gse_trx_all_5mo$p.adjust[37]

gseaplot2(Gse_trx2wk_B, geneSetID = 79, title = Gse_trx2wk_B$Description[79])
Gse_trx2wk_B$NES[79]
gseaplot2(Gse_trx, geneSetID = 49, title = Gse_trx$Description[49])
Gse_trx$NES[49]

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

row.names(Cnt_trx2wk_all) <- gsub(pattern = "[.].*", replacement = "", x= row.names(Cnt_trx2wk_all))
DE_trx2wk_DDR <- merge(Genelist_pDNAdamage_m,Cnt_trx2wk_all, by=0)
row.names(DE_trx2wk_DDR) <- DE_trx2wk_DDR[,4]
DE_trx2wk_DDR <- DE_trx2wk_DDR[,-1:-2]
DE_trx2wk_DDR$sumwt <- rowSums(DE_trx2wk_DDR[,7:8])
DE_trx2wk_DDR$sumhet <- rowSums(DE_trx2wk_DDR[,9:10])
DE_trx2wk_DDR_f <- filter(DE_trx2wk_DDR, sumhet - sumwt >500)
DE_trx2wk_DDR_f <- filter(DE_trx2wk_DDR_f, DE_trx2wk_DDR_f[,9] - DE_trx2wk_DDR_f[,8] > 0)
which(row.names(DE_trx2wk_DDR_f) =='Fus')
DE_trx2wk_DDR_f <- DE_trx2wk_DDR_f[-9,]

pheatmap(log2(DE_trx2wk_DDR_f[,7:10] + 1), cluster_cols = F, cellwidth = 40, cellheight = 14,border_color = NA,treeheight_row = 15, fontsize_col = 2,fontsize_row = 14, legend_breaks = c(-1.5,0,1.5), scale = 'row')
pheatmap(log2(DE_trx2wk_DDR_f[,12:15] + 1), cluster_cols = F, cellwidth = 40, cellheight = 14,border_color = NA,treeheight_row = 15, fontsize_col = 2,fontsize_row = 14, legend_breaks = c(-1.5,0,1.5), scale = 'row')

row<-c(8,12,4,14,6,9,7,3,5,10,19,11,2,4,18,1,15,16,17)
pheatmap(log2(DE_trx2wk_DDR_f[row,7:10] + 1), cluster_cols = F,cluster_rows = F, cellwidth = 40, cellheight = 14,border_color = NA,treeheight_row = 15, fontsize_col = 2,fontsize_row = 14, legend_breaks = c(-1.5,0,1.5), scale = 'row')
pheatmap(log2(DE_trx2wk_DDR_f[row,12:15] + 1), cluster_cols = F, cluster_rows = F, cellwidth = 40, cellheight =14,border_color = NA,treeheight_row = 15, fontsize_col = 5,fontsize_row = 14, legend_breaks = c(-1.5,0,1.5), scale = 'row')



write.csv(DE_doxo_DDR, "~/Desktop/DE_DDR_doxo.csv")





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

