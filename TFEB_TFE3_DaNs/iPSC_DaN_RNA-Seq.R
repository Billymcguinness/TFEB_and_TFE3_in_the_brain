library(BiocManager)
library(rnaseqGene)
library(GO.db)
library(clusterProfiler)
library(GenomeInfoDbData)
samples <- read.table(file.path("./DaN_Samples_meta.txt"), header =T)
samples
dir <- "./salmon/quants"
files <- file.path(dir, samples$runname, "quant.sf")
names(files) <- samples$run
files
all(file.exists(files))
samples
txdb<- GenomicFeatures::makeTxDbFromGFF("./Homo_sapiens.GRCh38.109.chr.gtf")
columns(txdb)
k <- keys(txdb, keytype ="TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID" ,keytype = "TXNAME")
head(tx2gene)
library(tximport)
library(tximeta)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = T)
Raw_counts <- txi$counts
write.table(Raw_counts, file = 'Raw_counts All samples.txt', sep = "\t", col.names = T, row.names = T)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ Cell_Line + Vector)
dim(ddsTxi)
head(row.names(ddsTxi))
ddsTxi$Vector <- relevel(ddsTxi$Vector , ref = "mApple")
levels(ddsTxi$Vector)
ddseq <- DESeq(ddsTxi)
fpkmDESeq2 <- fpkm(ddseq, robust = T)
head(fpkmDESeq2)
write.table(fpkmDESeq2, file = "All_samples_FPKM.txt", sep = "\t", row.names = T, col.names = T)
nrow(ddseq)
keepddseq <- rowSums(counts(ddseq)) > 30
ddseq <- ddseq[keepddseq,]
nrow(ddseq)
#16655 genes remaining once removed counts < 30
#TFE3 overexpression results
resT3mA <- results(ddseq, contrast = c("Vector", "TFE3", "mApple"))
head(resT3mA)
resT3mA <- resT3mA[order(resT3mA$pvalue),]
resT3mA
mcols(resT3mA, use.names = T)
summary(resT3mA)
#1041 up, 621 down
write.csv(resT3mA, file = "TFE3_DESeq2_output.csv", row.names = T, col.names = T)
TFE3_DaN_sigs <- na.omit(resT3mA)
TFE3_DaN_sigs <- TFE3_DaN_sigs[TFE3_DaN_sigs$padj < 0.05,]
TFE3_DaN_sigs[TFE3_DaN_sigs$padj < 0.05,]
summary(TFE3_DaN_sigs)
#805 genes up, 399 genes down
#TFEB overexpression results
resTBmA <- results(ddseq, contrast = c("Vector", "TFEB", "mApple"))
head(resTBmA)
resTBmA <- resTBmA[order(resTBmA$pvalue),]
resTBmA
mcols(resTBmA, use.names = T)
summary(resTBmA)
### 2303 upregulated genes, 2389 downregulated genes...
TFEB_DaN_sigs <- na.omit(resTBmA)
TFEB_DaN_sigs <- TFEB_DaN_sigs[TFEB_DaN_sigs$padj < 0.05,]
TFEB_DaN_sigs[TFEB_DaN_sigs$padj < 0.05,]
summary(TFEB_DaN_sigs)
###1965 upregulated genes, 1944 downregulated genes
## PCA plots
rld <- rlog(ddseq, blind = F)
head(assay(rld), 3)
plotPCA(rld, intgroup = c("Name", "Vector")) + geom_text(aes(label=c(samples$Name)), vjust=2) + geom_text(aes(label=c(samples$Vector)), vjust=-1)
##PC1 (43% variance) separates lines, PC2 (33% variance) separates by vector most effectively. 
## Gene ontology analysis
###BiocManager::install('AnnotationDbi', force = T)
###BiocManager::install('Biobase', force = T)
library(clusterProfiler)
T3 <- rownames(TFE3_DaN_sigs[TFE3_DaN_sigs$log2FoldChange > 0.4, ])
GO_results_TFE3 <- enrichGO(gene = T3, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
GO_results_TFE3
dotplot(GO_results_TFE3, showCategory = 20) + ggtitle("Top Biological Process GO - TFE3 overexpression")
### TFEB
TB <- rownames(TFEB_DaN_sigs[TFEB_DaN_sigs$log2FoldChange > 0.4, ])
GO_results_TFEB <- enrichGO(gene = TB, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
GO_results_TFEB
dotplot(GO_results_TFEB, showCategory = 20) + ggtitle("Top Biological Process GO - TFEB overexpression")