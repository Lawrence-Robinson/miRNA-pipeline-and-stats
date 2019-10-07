#========================# STEPHAN'S SECTION #============================#
biocLite("TCC")
library(TCC); library(DESeq); library(reshape2); library(data.table); library(ggplot2); library(ggrepel); library(viridis); library(baySeq); library(DESeq2); library(pasilla);
#setwd("G:/rp2/")
setwd("D:\\Users\\Lawrence miRNA")

# Load data
data = read.csv("mirdeep_out.csv")
summary = read.csv("resultX.csv", sep = "\t")

# Isolate mirna, rfam and ranfold columns from summary
summary = summary[substring(as.vector(summary[,"miRDeep2.score"]), 1, nchar("chr")) == "chr", c("miRDeep2.score", "novel.miRNAs..estimated.true.positives", "excision.gearing")]
colnames(summary) = c("id", "rfam_alert", "significant_randfold_p-value")

# Isolate a unique 'id' list of where novel = 0.
mir_mirbase = unique(data[data$mirbase_id != "na",c("id","mirbase_id", "mirbase_accession")])

# Remove these novel = 0 identifiers from csv file.
mir_count_all = data[data$mirbase_id == "na",]

# Fill the missing column with the known miRNA mirbase ID
mir_count_all$mirbase_id = apply(mir_count_all, 1, function(x) ifelse(any(mir_mirbase$id %in% x[1]), as.character(mir_mirbase[which(mir_mirbase$id %in% x[1]),"mirbase_id"]), "none"))
mir_count_all$mirbase_accession = apply(mir_count_all, 1, function(x) ifelse(any(mir_mirbase$id %in% x[1]), as.character(mir_mirbase[which(mir_mirbase$id %in% x[1]),"mirbase_accession"]), "none"))
#LB: BTW if you have 1:many mapping, this ignores them / collates them into the alphabetically first ID - is this intentional?

# Aggregate the data set to get read totals across sample
mir_count_sum = aggregate(count ~ id + sample + group + region + score + exp_seq + pri_seq + mirbase_id + mirbase_accession + related_accession + related_id, 
                          data = mir_count_all, 
                          FUN = "sum")

# Add rfam alert and randfold p value significance
mir_count_sum$rfam_alert = sapply(as.vector(mir_count_sum$id), FUN = function(x) ifelse(x %in% summary[summary$rfam_alert != "-", "id"], as.character(summary[summary$id == x, "rfam_alert"]), "-"))
mir_count_sum$sig._randfold_p = sapply(as.vector(mir_count_sum$id), FUN = function(x) ifelse(x %in% summary[summary$`significant_randfold_p-value` == "yes", "id"], "yes", "no"))

# Generate unique id
mir_count_sum$unique_id = paste(mir_count_sum$id, mir_count_sum$sample, sep = "-")

#isolate only reads that are mature, have a score of greater than 4 and no rfam alerts
mir_count_mat = mir_count_sum[mir_count_sum$score >= 4 & 
                                mir_count_sum$region == "mature" &
                                mir_count_sum$rfam_alert == "-",]

#===================== LAWRENCE SECTION =======================# (rest of file) #

mir_count_unknown = mir_count_mat_sub[,-(2:5)]
mir_count_unknown = mir_count_unknown[,-(8:11)]
mir_count_unknown = mir_count_unknown[,-(2:3)]
mir_count_unknown = mir_count_unknown[!duplicated(mir_count_unknown),]
mir_count_unknown = mir_count_unknown[mir_count_unknown$mirbase_accession == "none" & mir_count_unknown$related_accession == "none",]

mir_count_week =mir_count_mat_sub[, -(2)] 
mir_count_week =mir_count_week[, -(11:14)]
mir_count_week =mir_count_week[, -(4:6)]
mir_count_week =mir_count_week[mir_count_week$group =="One_Week",]
mir_count_week = mir_count_week[!duplicated(mir_count_week),]
mir_count_week = mir_count_week[mir_count_week$mirbase_accession == "none" & mir_count_week$related_accession == "none",]

mir_count_month =mir_count_mat_sub[, -(2)] 
mir_count_month =mir_count_month[, -(11:14)]
mir_count_month =mir_count_month[, -(4:6)]
mir_count_month =mir_count_month[mir_count_month$group =="One_Month",]
mir_count_month = mir_count_month[!duplicated(mir_count_month),]
mir_count_month = mir_count_month[mir_count_month$mirbase_accession == "none" & mir_count_month$related_accession == "none",]

mir_count_three =mir_count_mat_sub[, -(2)] 
mir_count_three =mir_count_three[, -(11:14)]
mir_count_three =mir_count_three[, -(4:6)]
mir_count_three =mir_count_three[mir_count_three$group =="Three_Month",]
mir_count_three = mir_count_three[!duplicated(mir_count_three),]
mir_count_three = mir_count_three[mir_count_three$mirbase_accession == "none" & mir_count_three$related_accession == "none",]

mir_count_adult =mir_count_mat_sub[, -(2)] 
mir_count_adult =mir_count_adult[, -(11:14)]
mir_count_adult =mir_count_adult[, -(4:6)]
mir_count_adult =mir_count_adult[mir_count_adult$group =="Adult",]
mir_count_adult = mir_count_adult[!duplicated(mir_count_adult),]
mir_count_adult = mir_count_adult[mir_count_adult$mirbase_accession == "none" & mir_count_adult$related_accession == "none",]

count_all = merge(mir_count_adult, mir_count_month, all=TRUE)
count_all = merge(mir_count_week, mir_count_three, all=TRUE)
count_uniq = count_all[, -(2)]

count_uniq = count_uniq[!duplicated(count_uniq),]

rm(data, mir_count_all, mir_count_sum, mir_mirbase, summary)

### Run differential expression analysis on different subsets (pairwise etc.) ###
norm_counts = list()
all_results = NULL

mir_count_mat_sub = mir_count_mat
count_table = acast(mir_count_mat_sub, id ~ sample, value.var = "count", fill = 0)
count_table = count_table[,order(colnames(count_table), decreasing = TRUE)]
print(nrow(count_table))
#no need to normalise the counts as DESeq will do this for us
#we want to aggregate the counts per group for differential expression after normalisation 

library(apeglm); library(vsn); library(pheatmap); library(RColorBrewer); library(viridis); library(viridisLite); library(ggplot2)

sample = substr(colnames(count_table),1,3)
a = "XAdult"
m = "One_Month"
t = "Three_Month"
w = "A_Week"
group=c(w,w,w,w,w,w,w,t,t,t,t,t,t,t,t,m,m,m,m,m,m,m,a,a,a,a,a,a,a,a)
#group = c(a,a,a,a,a,a,a,a,m,m,m,m,m,m,m,m,t,t,t,t,t,t,t,w,w,w,w,w,w,w)
coldata= cbind(sample, group)
count_table2 = count_table[,order(colnames(count_table), decreasing = TRUE)]
coldata2 = coldata
cts = count_table2

######################################################## performa a DESeq with just W and M ####################################################################

ct1 = count_table2[,-(23:30)]
count_tableWM = ct1[,-(8:15)]
groupWM = c(w,w,w,w,w,w,w,m,m,m,m,m,m,m)
sampleWM =substr(colnames(count_tableWM),1,3) 
coldataWM = cbind(sampleWM, groupWM)

ddsWM <- DESeqDataSetFromMatrix(countData = count_tableWM, colData = coldataWM, design = ~ groupWM)
keep <- rowSums(counts(ddsWM)) >=100
#ddsWM <- ddsWM[keep,]
ddsWM <- DESeq(ddsWM)
res_WM <- results(ddsWM, contrast=c("groupWM", "One_Month", "A_Week"))
resultsNames(ddsWM)

####################################################### DESeq with One Week and Three Month #################################################################

count_tableWT = count_table2[,-(16:30)]
groupWT = c(w,w,w,w,w,w,w,t,t,t,t,t,t,t,t)
sampleWT = substr(colnames(count_tableWT),1,3)
coldataWT = cbind(sampleWT, groupWT)

ddsWT <- DESeqDataSetFromMatrix(countData = count_tableWT, colData = coldataWT, design = ~ groupWT)
keep <- rowSums(counts(ddsWT)) >=100
#ddsWT <- ddsWT[keep,]
ddsWT <- DESeq(ddsWT)
res_WT <- results(ddsWT, contrast=c("groupWT", "Three_Month", "A_Week"))
resultsNames(ddsWT)

####################################################### DESeq with One Week and Adult #########################################################################

count_tableWA = count_table2[,-(8:22)]
groupWA = c(w,w,w,w,w,w,w,a,a,a,a,a,a,a,a)
sampleWA = substr(colnames(count_tableWA),1,3)
coldataWA = cbind(sampleWA, groupWA)

ddsWA <- DESeqDataSetFromMatrix(countData = count_tableWA, colData = coldataWA, design = ~ groupWA)
#keep <- rowSums(counts(ddsWA)) >= 50
#ddsWA <- ddsWA[keep,]
ddsWA <- DESeq(ddsWA)
res_WA <- results(ddsWA, contrast=c("groupWA", "XAdult","A_Week"))
resultsNames(ddsWA)

################################################### All Results Stats handling ################################################################################

resOrdered_WM <- res_WM[order(res_WM$padj),]
resOrdered_WT <- res_WT[order(res_WT$padj),]
resOrdered_WA <- res_WA[order(res_WA$padj),]
write.csv(as.data.frame(resOrdered_WA), file="DESeqWk_vs_Adult.csv")
write.csv(as.data.frame(resOrdered_WT), file="DESeqWk_vs_Three.csv")
write.csv(as.data.frame(resOrdered_WM), file="DESeqWk_vs_Month.csv")
summary(res_WM)
summary(res_WT)
summary(res_WA)

hist(res_WA$padj, breaks=20, col="grey" )
hist(res_WT$padj, breaks=20, col="grey" )
hist(res_WM$padj, breaks=20, col="grey" )

sum(res_WA$log2FoldChange > 2, na.rm = TRUE)
sum(res_WA$log2FoldChange < -2, na.rm = TRUE)

sum(res_WT$log2FoldChange > 0 & res_WA$padj < 0.05, na.rm = TRUE)
sum(res_WT$log2FoldChange < 0 & res_WA$padj < 0.05, na.rm = TRUE)

sum(res_WM$log2FoldChange > 0, na.rm = TRUE)
sum(res_WM$log2FoldChange < 0 & res_WA$padj < 0.05, na.rm = TRUE)

sum(res_WA$padj < 0.05, na.rm=TRUE)
sum(res_WT$padj < 0.05, na.rm=TRUE)
sum(res_WM$padj < 0.05, na.rm=TRUE)

plotMA(res_WM, ylim=c(-2, 2))
plotMA(res_WT, ylim=c(-2, 2))
plotMA(res_WA, ylim=c(-2, 2))

res_WALFC <- lfcShrink(ddsWA, coef="groupWA_XAdult_vs_A_Week", type="apeglm")
res_WTLFC <- lfcShrink(ddsWT, coef="groupWT_Three_Month_vs_A_Week", type="apeglm")
res_WMLFC <- lfcShrink(ddsWM, coef="groupWM_One_Month_vs_A_Week", type="apeglm")

resApeT <- lfcShrink(ddsWM, coef=2, type="apeglm", lfcThreshold=1)
plotMA(resApeT, ylim=c(-5,5), cex=1)
abline(h=c(-1,1), col="dodgerblue", lwd=2)

resApeT <- lfcShrink(ddsWT, coef=2, type="apeglm", lfcThreshold=1)
plotMA(resApeT, ylim=c(-5,5), cex=1)
abline(h=c(-1,1), col="dodgerblue", lwd=2)

resApeT <- lfcShrink(ddsWA, coef=2, type="apeglm", lfcThreshold=1)
plotMA(resApeT, ylim=c(-5,5), cex=1)
abline(h=c(-1,1), col="dodgerblue", lwd=2)

# independant filtering
metadata(res_WA)$filterThreshold

plot(metadata(res_WM)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res_WM)$lo.fit, col="red")
abline(v=metadata(res_WM)$filterTheta, col="blue")

par(mar=c(8,5,2,2))
boxplot(log10(assays(ddsWA)[["cooks"]]), range=0, las=2)
par(mar=c(8,5,2,2))
boxplot(log10(assays(ddsWT)[["cooks"]]), range=0, las=2)
par(mar=c(8,5,2,2))
boxplot(log10(assays(ddsWM)[["cooks"]]), range=0, las=2)

plotMA(res_WALFC, ylim=c(-2,2))
plotMA(res_WTLFC, ylim=c(-2,2))
plotMA(res_WMLFC, ylim=c(-2,2))

# graphical results Week v Adult
breaksList = seq(0, 40, by = 1)
vsd <- varianceStabilizingTransformation(ddsWA, blind=FALSE)
rld <- rlog(ddsWA, blind=FALSE)

plotCounts(ddsWA, gene=which.min(res_WA$padj), intgroup="groupWA")
ntd <- normTransform(ddsWA)
meanSdPlot(assay(ntd))
select <- order(rowMeans(counts(ddsWA,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddsWA)[,c("groupWA")])
pheatmap(assay(ntd[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df))
#heatmap
sampleDists <- dist(t(assay(ntd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(ntd$sampleWA, sep="-" )
colnames(sampleDistMatrix) <- paste(ntd$sampleWA, sep="-" )
colors <- colorRampPalette(brewer.pal(8, "Blues") )(150)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col=plasma(40), breaks = breaksList)

pcaData <- plotPCA(vsd, intgroup=c("groupWA", "sampleWA"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=groupWA, shape=groupWA)) + geom_point(size=5) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed()

plotDispEsts(ddsWA)

# graphical results Week v three month
vsd <- varianceStabilizingTransformation(ddsWT, blind=FALSE)
rld <- rlog(ddsWT, blind=FALSE)

plotCounts(ddsWT, gene=which.min(res_WT$padj), intgroup="groupWT")
ntd <- normTransform(ddsWT)
meanSdPlot(assay(ntd))
select <- order(rowMeans(counts(ddsWT,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddsWT)[,c("groupWT")])
pheatmap(assay(ntd[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df))
#heatmap
sampleDists <- dist(t(assay(ntd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(ntd$sampleWT, sep="-" )
colnames(sampleDistMatrix) <- paste(ntd$sampleWT, sep="-" )
colors <- colorRampPalette(brewer.pal(7, "Blues")) (255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col=plasma(40), breaks = breaksList)

pcaData <- plotPCA(vsd, intgroup=c("groupWT", "sampleWT"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=groupWT, shape = groupWT)) + geom_point(size=5) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed()

plotDispEsts(ddsWT)

#graphical results week v one month
vsdWM <- varianceStabilizingTransformation(ddsWM, blind=FALSE)
rldWM <- rlog(ddsWM, blind=FALSE)

plotCounts(ddsWM, gene=which.min(res_WM$padj), intgroup="groupWM")
ntd <- normTransform(ddsWM)
meanSdPlot(assay(ntd))
select <- order(rowMeans(counts(ddsWM,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(ddsWM)[,c("groupWM")])
pheatmap(assay(ntd[select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df))
#heatmap
sampleDists <- dist(t(assay(ntd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(ntd$sampleWM, sep="-" )
colnames(sampleDistMatrix) <- paste(ntd$sampleWM, sep="-" )
color <- viridis
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col=plasma(40), breaks = breaksList)

pcaData <- plotPCA(vsdWM, intgroup=c("groupWM", "sampleWM"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=groupWM, shape=groupWM)) + geom_point(size=5) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed()

plotDispEsts(ddsWM)
#### simple comparison of all the data points

ntd <- normTransform(dds)
sampleDists <- dist(t(assay(ntd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(ntd$sample, sep="-" )
colnames(sampleDistMatrix) <- paste(ntd$sample, sep="-" )
colors <- colorRampPalette(brewer.pal(7, "Blues")) (255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col=plasma(100))

vsd = varianceStabilizingTransformation(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group)) + geom_point(size=3) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) + coord_fixed()

plotDispEsts(dds)

######## add results to existing data frames ########

res_table = mir_count_mat_sub[mir_count_mat_sub$sig._randfold_p =="yes",]
res_table = res_table[,-(5:7)]
res_table = res_table[,-(9:11)]
res_table = res_table[,-(2)]
res_table = res_table[,-(3)]
res_WA.dt = read.csv("DESeqWk_vs_Adult.csv")
colnames(res_WA.dt)[1] <- "id"
RESULTS_WA = merge(res_table, res_WA.dt[, c("id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")], by = "id")
RESULTS_WA = RESULTS_WA[RESULTS_WA$padj <= 0.05 & RESULTS_WA$group == "One_Week",]
RESULTS_WA = na.omit(RESULTS_WA)
RESULTS_WA = RESULTS_WA[!duplicated(RESULTS_WA$id),]
RESULTS_WA_short1 = RESULTS_WA[RESULTS_WA$log2FoldChange <= 0 ,]
RESULTS_WA_short1 = RESULTS_WA[RESULTS_WA$mirbase_id != "none" & RESULTS_WA$related_id == "none",]
RESULTS_WA_short1 = RESULTS_WA[RESULTS_WA$mirbase_id == "none" & RESULTS_WA$related_id == "none",]
RESULTS_WA_short1 = RESULTS_WA[RESULTS_WA$related_id != "none",]


RESULTS_WA_short2 = RESULTS_WA[RESULTS_WA$log2FoldChange <= 0 ,]
RESULTS_WA_short = merge(RESULTS_WA_short1, RESULTS_WA_short2, all= TRUE)
RES_WA_TABLE <- RESULTS_WA_short[order(RESULTS_WA_short$log2FoldChange),]
RES_WA_TABLE = RES_WA_TABLE[,-(2:4)]
RES_WA_TABLE = RES_WA_TABLE[,-(4:5)]
RES_WA_TABLE = RES_WA_TABLE[,-(5:7)]
formattable(RES_WA_TABLE)


RESULTS_WA_miRNA = RESULTS_WA[RESULTS_WA$mirbase_id != "none" | RESULTS_WA$related_id != "none",]
humanWA = RESULTS_WA[RESULTS_WA$related_id != "none" & RESULTS_WA$mirbase_id == "none",]
sheepWA = RESULTS_WA[RESULTS_WA$mirbase_id != "none" & RESULTS_WA$related_id == "none",]
IPA_WA = merge(humanWA, sheepWA, all=TRUE)
IPA_WA = IPA_WA[, -(1:2)]
IPA_WA = IPA_WA[, -(5)]
write.csv(as.data.frame(IPA_WA), file = "IPA_WA.xls")

# week and three month
res_WT.dt = read.csv("DESeqWk_vs_Three.csv")
colnames(res_WT.dt)[1] <- "id"
RESULTS_WT = merge(res_table, res_WT.dt[, c("id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")], by = "id")
RESULTS_WT = RESULTS_WT[RESULTS_WT$padj <= 0.05 & RESULTS_WT$group == "One_Week",]
RESULTS_WT = na.omit(RESULTS_WT)
RESULTS_WT = RESULTS_WT[!duplicated(RESULTS_WT$id),]
RESULTS_WT_short1 = RESULTS_WT[RESULTS_WT$log2FoldChange < 0 ,]
RESULTS_WT_short1 = RESULTS_WT[RESULTS_WT$mirbase_id != "none" & RESULTS_WT$related_id == "none",]
RESULTS_WT_short1 = RESULTS_WT[RESULTS_WT$mirbase_id == "none" & RESULTS_WT$related_id == "none",]
RESULTS_WT_short1 = RESULTS_WT[RESULTS_WT$related_id != "none",]


RESULTS_WT_miRNA = RESULTS_WT[RESULTS_WT$mirbase_id != "none" | RESULTS_WT$related_id != "none",]
humanWT = RESULTS_WT[RESULTS_WT$related_id != "none" & RESULTS_WT$mirbase_id == "none",]
sheepWT = RESULTS_WT[RESULTS_WT$mirbase_id != "none" & RESULTS_WT$related_id == "none",]
IPA_WT = merge(humanWT, sheepWT, all=TRUE)
IPA_WT = IPA_WT[, -(1:2)]
IPA_WT = IPA_WT[, -(5)]
write.csv(as.data.frame(IPA_WT), file = "IPA_WT.xls")

#week and one month
res_WM.dt = read.csv("DESeqWk_vs_Month.csv")
colnames(res_WM.dt)[1] <- "id"
RESULTS_WM = merge(res_table, res_WM.dt[, c("id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")], by = "id")
RESULTS_WM = RESULTS_WM[RESULTS_WM$padj <= 0.05 & RESULTS_WM$group == "One_Week",]
RESULTS_WM = na.omit(RESULTS_WM)
RESULTS_WM = RESULTS_WM[!duplicated(RESULTS_WM$id),]
RESULTS_WM_short1 = RESULTS_WM[RESULTS_WM$log2FoldChange > 0 ,]
RESULTS_WM_short1 = RESULTS_WM[RESULTS_WM$mirbase_id != "none" & RESULTS_WM$related_id == "none",]
RESULTS_WM_short1 = RESULTS_WM[RESULTS_WM$mirbase_id == "none" & RESULTS_WM$related_id == "none",]
RESULTS_WM_short1 = RESULTS_WM[RESULTS_WM$related_id != "none",]


RESULTS_WM_miRNA = RESULTS_WM[RESULTS_WM$mirbase_id != "none" | RESULTS_WM$related_id != "none",]
humanWM = RESULTS_WM[RESULTS_WM$related_id != "none" & RESULTS_WM$mirbase_id == "none",]
sheepWM = RESULTS_WM[RESULTS_WM$mirbase_id != "none" & RESULTS_WM$related_id == "none",]
IPA_WM = merge(humanWM, sheepWM, all=TRUE)
IPA_WM = IPA_WM[, -(1:2)]
IPA_WM = IPA_WM[, -(5)]
write.csv(as.data.frame(IPA_WM), file = "IPA_WM.xls")

###### post IPA rearrangements ########

library("formattable"); library("dplyr"); library("tidyr"); library(knitr)
IPA_OUT_WA = read.csv("WA_IPA_OUT.csv")
IPA_OUT_WA = aggregate(
  IPA_OUT_WA$Symbol.1 ~ IPA_OUT_WA$ID + IPA_OUT_WA$Symbol + IPA_OUT_WA$Expr_Log_Ratio, data=IPA_OUT_WA, FUN=paste, collapse=' '
)
names(IPA_OUT_WA)[1]<-"Accession"
names(IPA_OUT_WA)[2]<-"miRbase ID"
names(IPA_OUT_WA)[3]<-"Log2Fold Change"
names(IPA_OUT_WA)[4]<-"Gene Targets"
WA_TARGETS_G1 = IPA_OUT_WA[IPA_OUT_WA$`Log2Fold Change` >= 1,]
WA_TARGETS_L1 = IPA_OUT_WA[IPA_OUT_WA$`Log2Fold Change`<= -1,]
WA_TARGETS = merge(WA_TARGETS_G1, WA_TARGETS_L1, all=TRUE)
WA_TARGETS[sort(WA_TARGETS$`Log2Fold Change`),]
WA_TARGETS
sort() help

formattable(WA_TARGETS)
write.table((WA_TARGETS), file="WA_TARGETS.txt", quote = FALSE, row.names = F)


###### week v three month #######











