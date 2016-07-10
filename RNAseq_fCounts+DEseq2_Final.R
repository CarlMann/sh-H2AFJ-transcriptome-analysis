###########################################################################
###      Install Bioconductor and required packages if necessary   ########
##########################################################################
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("Rsubread", "DESeq2", "RColorBrewer", "pheatmap", "AnnotationDbi", "org.Hs.eg.db", "genefilter" )) 

# Download ArrayExpress files (accession number E-MTAB-4920) containing fastq paired-end DNA sequences. 
# Use trim-galore(with Cutadapt) and FastQC to remove Illumina TruSeq adapters
# and control quality. Use Tophat2 to align sequences to the hg19 genome.
# Transfer the tophat2 accepted_hits.bam file for the 3 reps of sh-NT and sh3-H2AFJ to 
# a sub-directory labelled "DataBAM" in the R working directory.

# Use featureCounts and inbuilt hg19 annotation table to build count table:
# paired-end reads specified, so featureCount will only align the paired-seqs 
# as a single fragment. Default option summarizes counts to MetaFeatures (genes composed
# of combined exons, rather than counting to individual exons). Default allowMultiOverlap=FALSE
# means that only uniquely aligned reads will be counted.

library("Rsubread")

# accepted_hits.bam files of paired-end sequences mapped to the hg19 genome, produced by Tophat2,
# should be in a DataBam sub-directory within the working directory
bamFiles <- list.files(path = file.path(getwd(),"DataBAM"), pattern = ".bam", full.names=TRUE)
fc <- featureCounts(files=c(bamFiles), nthreads=2,isPairedEnd=TRUE,annot.inbuilt="hg19")

# Feature Counts to DESeq2
# DEseq2 analysis following Mike Love: http://master.bioconductor.org/help/workflows/rnaseqGene/
fc2 <- fc$counts
colnames(fc2) <- c("NTSen4","NTSen5","NTSen6","shSEN4","shSEN5","shSEN6")
Sample <- c("NTSen4","NTSen5","NTSen6","shSEN4","shSEN5","shSEN6")
Condition <- c("NT","NT","NT","shJ","shJ","shJ")
(sampleTable <- cbind(Sample, Condition)) #create and display
coldata <- as.data.frame(sampleTable)
library("DESeq2")
(ddsMat <- DESeqDataSetFromMatrix(countData = fc2, colData = coldata, design = ~ Condition)) #create and display summary

round( colSums(assay(ddsMat)) / 1e6, 1 )  # number of paired-end  sequence fragments per million mapped to genes for each sample:
#      1    2    3    4    5    6 
#     17.2 16.5  9.5 20.2 14.4 14.4 

# Exploratory analysis and visualization for DESeq2 analysis:
# 2 types of anal: 
# -transformations of the counts in order to visually explore sample relationships
# -back to the original raw counts for statistical testing
# This is critical because the statistical testing methods rely on original count data 
# (not scaled or transformed) for calculating the precision of measurements.

# Filter out genes (rows) that have no counts, or only a single count across all samples:

dds2 <- ddsMat[rowSums(counts(ddsMat)) > 1,]
nrow(dds2)  # 20,083 genes with at least 2 counts in at least one sample

#Many common statistical methods for exploratory analysis of multidimensional data, for example 
# clustering and principal components analysis (PCA), work best for data that generally has the 
# same range of variance at different ranges of the mean values. For RNA-seq raw counts, however, 
# the variance grows with the mean. Use regularized-logarithm transformation or rlog 
# (Love, Huber, and Anders 2014) to stabilize the variance across the mean. 

rld2 <- rlog(dds2, blind=FALSE)

# Use the R function dist to calculate the Euclidean distance between samples:
sampleDists2 <- dist( t( assay(rld2) ) )
sampleDists2
library("pheatmap")
library("RColorBrewer")

# Provide sampleDists to the clustering_distance argument of the pheatmap function. 
# Otherwise the pheatmap function would assume that the matrix contains the data values themselves, 
# and would calculate distances between the rows/columns of the distance matrix, which is not desired.

sampleDistMatrix2 <- as.matrix( sampleDists2 )
rownames(sampleDistMatrix2) <- rld2$Sample
colnames(sampleDistMatrix2) <- rld2$Sample
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix2,
         clustering_distance_rows=sampleDists2,
         clustering_distance_cols=sampleDists2,
         col=colors)
# NTSen4 appears as an outlier

# Differential expression pipeline after featureCounts:
dds2 <- DESeq(dds2)
(res2 <- results(dds2)) # create and display
summary(res2)
resSig2 <- subset(res2, padj < 0.1)

library("AnnotationDbi")
library("org.Hs.eg.db")

resSig2$symbol <- mapIds(org.Hs.eg.db,keys=row.names(resSig2),column="SYMBOL",keytype="ENTREZID",multiVals="first")
resSig2 <- resSig2[ order(resSig2$log2FoldChange), ]

library("genefilter")
siggenes2 <- rownames(resSig2)
mat2 <- assay(rld2)[siggenes2,]
mat2<- mat2- rowMeans(mat2)
namevec2 <- mapIds(org.Hs.eg.db,keys=row.names(mat2),column="SYMBOL",keytype="ENTREZID",multiVals="first")
rownames(mat2) <- mapIds(org.Hs.eg.db,keys=row.names(mat2),column="SYMBOL",keytype="ENTREZID",multiVals="first")
df2 <- as.data.frame(colData(rld2)[,c("Sample","Condition")])
pheatmap(mat2, annotation_col=df2, cluster_rows=F,cellheight = 10) # complete heat map 
# using differentially expressed genes, we observe the expected 2 distinct clusters: sh2-H2AFJ and sh-NoTarget
RNAseqDE <- as.data.frame(resSig2)
write.csv(RNAseqDE, file = "RNAseqDE_fCounts+DESeq2.csv")

# produce condensed heat map of genes significantly down-regulated by sh2-H2AFJ knockdown versus sh-NT in senescence
DownRNAseq <- RNAseqDE[which(RNAseqDE$log2FoldChange<0),]
RNAseqMat <- as.matrix(DownRNAseq$log2FoldChange)
rownames(RNAseqMat) <- DownRNAseq$symbol
pheatmap(RNAseqMat,scale='none',cluster_cols=F, cluster_rows=F,show_rownames = T,cellwidth = 30,cellheight = 19.5,
         col=brewer.pal(9,"YlOrRd"))
         
# FKPM Table
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")]) #create DGEList with Rsubread
dim(x$counts)        	# 25,702 genes
fkpmTable <- rpkm(fc$counts, fc$annotation$Length)
fkpmTable <- fkpmTable[rowSums(fkpmTable)>1,] #14,574 genes with sum greater than 1 FKPM for all samples
fkpmTable <- as.data.frame(fkpmTable) # convert to data frame so that gene symbol column can be added
symbol <- mapIds(org.Hs.eg.db,keys=row.names(fkpmTable),column="SYMBOL",keytype="ENTREZID",multiVals="first") # get gene symbols 
symbol <- as.vector(symbol) # convert to vector of gene symbols
fkpmTable$Symbol <- symbol # create column of gene symbols
fkpmTable <- fkpmTable[order(-fkpmTable$NTSen5),] # order by decreasing FKPM for NTSen5 sample
write.csv(fkpmTable, file="sh2_shNT_SEN_RNAseq_FKPM.csv") # write table for publication

# sessionInfo()
# R version 3.2.2 (2015-08-14)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
# Running under: OS X 10.11.5 (El Capitan)
# 
# locale:
#         [1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8
# 
# attached base packages:
#         [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#         [1] pheatmap_1.0.8             genefilter_1.52.1          org.Hs.eg.db_3.2.3         RSQLite_1.0.0             
# [5] DBI_0.4-1                  AnnotationDbi_1.32.3       DESeq2_1.10.1              RcppArmadillo_0.7.100.3.1 
# [9] Rcpp_0.12.5                SummarizedExperiment_1.0.2 GenomicRanges_1.22.4       GenomeInfoDb_1.6.3        
# [13] IRanges_2.4.8              S4Vectors_0.8.11           Biobase_2.30.0             BiocGenerics_0.16.1       
# [17] RColorBrewer_1.1-2         BiocInstaller_1.20.3       Rsubread_1.20.6            edgeR_3.12.1              
# [21] limma_3.26.9              
# 
# loaded via a namespace (and not attached):
#         [1] futile.logger_1.4.1  plyr_1.8.4           XVector_0.10.0       futile.options_1.0.0 tools_3.2.2         
# [6] zlibbioc_1.16.0      rpart_4.1-10         annotate_1.48.0      gtable_0.2.0         lattice_0.20-33     
# [11] Matrix_1.2-6         gridExtra_2.2.1      cluster_2.0.4        locfit_1.5-9.1       nnet_7.3-12         
# [16] grid_3.2.2           data.table_1.9.6     XML_3.98-1.4         survival_2.39-5      BiocParallel_1.4.3  
# [21] foreign_0.8-66       latticeExtra_0.6-28  Formula_1.2-1        geneplotter_1.48.0   ggplot2_2.1.0       
# [26] lambda.r_1.1.7       Hmisc_3.17-4         scales_0.4.0         splines_3.2.2        xtable_1.8-2        
# [31] colorspace_1.2-6     acepack_1.3-3.3      munsell_0.4.3        chron_2.3-47     
