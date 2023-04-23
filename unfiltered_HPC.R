########################################################################
# [First-time use only] Install DESeq2 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2", dependencies=TRUE)


########################################################################
# Sec 1: data input; generate DESeqDataSet

# load DESeq2 (each time to use it)
library(DESeq2)

# the input count table "FinalProjectCounts.txt" should be in your current working directory (or you should provide the full path to it)
##readcounts <- read.csv("REF_featureCounts.csv")
##dim(readcounts) #check the table size
##head(readcounts, 10)
##View(readcounts)

unfilteredcounts <- read.table('counts_unfiltered.txt', sep="\t", header=TRUE)
unfilteredcounts$Geneid <- gsub("\\..*", "", unfilteredcounts$Geneid)
dim(unfilteredcounts)
head(unfilteredcounts, 10)
View(unfilteredcounts)


#Step 2: Generate the count table
# Combine the counts columns in finalcounts into a count table named final.count.table
# Set the column names as "HW1_count" and "demo_count", respectively 
# Set the row names to individual Geneid (note: Geneid is in column 1 of either count file)
# Check and show the first 10 rows of final.count.table

unfilteredcounts2 <- subset(unfilteredcounts, select=c(7:18))
View(unfilteredcounts2)

unfiltered.count.table <- data.frame(unfilteredcounts2)
View(unfiltered.count.table)
names(unfiltered.count.table)

#I couldn't figure out an easier way without error messages so... just run these all.

colnames(unfiltered.count.table)[colnames(unfiltered.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.unfiltered.SRR349741.fastq.gz_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.1"
colnames(unfiltered.count.table)[colnames(unfiltered.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.unfiltered.SRR349742.fastq.gz_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.2"
colnames(unfiltered.count.table)[colnames(unfiltered.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.unfiltered.SRR349743.fastq.gz_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.3"
colnames(unfiltered.count.table)[colnames(unfiltered.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.unfiltered.SRR349744.fastq.gz_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.4"
colnames(unfiltered.count.table)[colnames(unfiltered.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.unfiltered.SRR349745.fastq.gz_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.5"
colnames(unfiltered.count.table)[colnames(unfiltered.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.unfiltered.SRR349746.fastq.gz_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.6"
colnames(unfiltered.count.table)[colnames(unfiltered.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.unfiltered.SRR349747.fastq.gz_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.7"
colnames(unfiltered.count.table)[colnames(unfiltered.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.unfiltered.SRR349748.fastq.gz_Aligned.sortedByCoord.out.bam"] ="Non.tumor_sample.1"
colnames(unfiltered.count.table)[colnames(unfiltered.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.unfiltered.SRR349749.fastq.gz_Aligned.sortedByCoord.out.bam"] ="Non.tumor_sample.2"
colnames(unfiltered.count.table)[colnames(unfiltered.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.unfiltered.SRR349750.fastq.gz_Aligned.sortedByCoord.out.bam"] ="Non.tumor_sample.3"
colnames(unfiltered.count.table)[colnames(unfiltered.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.unfiltered.SRR349751.fastq.gz_Aligned.sortedByCoord.out.bam"] ="Non.tumor_sample.4"
colnames(unfiltered.count.table)[colnames(unfiltered.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.unfiltered.SRR349752.fastq.gz_Aligned.sortedByCoord.out.bam"] ="Non.tumor_sample.5"

names(unfiltered.count.table)
View(unfiltered.count.table)

### Input data cleanup 
# Get and keep a "numeric matrix" (to be iput into countData of the DESeqDataSet object)
# remove columns (1&2) that do not contain read counts (after saving the gene IDs $gene)

gene_id_unique <- unfilteredcounts$Geneid
gene_id_unique <- make.names(gene_id_unique, unique=TRUE)
row.names(unfiltered.count.table) <- gene_id_unique  #row.names(): gets/sets row names (vector) of a data frame
dim(unfiltered.count.table) #check the result after manipulating your data!
View(unfiltered.count.table) #note that row numbers are replaced with row names


### generate the DESeqDataSet object (DESeq2-specific R object)
# DESeqDataSet has 3 fields: countData, ColData, design

# make a data frame (colData, the meta-data) where row.names match the column names (sample names) in countData 
sample.group <- gsub("\\.[0-9]+", "", names(unfiltered.count.table))
View(sample.group)
#sample.group <-names(final.count.table)
unfiltered_info <- data.frame(row.names=names(unfiltered.count.table), condition=sample.group)
unfiltered_info #check the result
View(unfiltered_info)

# <function_name>(<argname>=<value>, .), "=" signals named argument passing in a function call
### "condition" sets a vector that stores sample info to be used in the "design" field in DESeqDataSet (see below)  
# generate the DESeqDataSet (with required fields: countData, ColData, design)

unfiltered.ds <- DESeqDataSetFromMatrix (countData = unfiltered.count.table, colData = unfiltered_info, design = ~condition) #design refers to a column(s) in colData (here, ~condition indicates the "condition" column) 


########################################################################
# Sec 2: rlog normalization/transformation & subsequent QC analysis/plots (correlation, clustering, PCA)

### Normalization using rlog
###
# calculate the size factor and add it to the data set
unfiltered.ds <- estimateSizeFactors (unfiltered.ds)
# if you check colData() again, you can see it now contains the sizeFactors
colData (unfiltered.ds)

# rlog() function returns values that are both normalized (by applying the size factor) and transformed to the log2 scale; it also shrinks the variance of low read counts
unfiltered_DESeq.rlog <- rlog(unfiltered.ds, blind = TRUE) # use blind = FALSE if large differences exist in a large proportion of genes.
unfiltered_rlog.counts <- assay(unfiltered_DESeq.rlog) # assay() function gets the read counts (integer) or normalized values (decimal), depending on the argument value
head(unfiltered_rlog.counts) # check and compare to head(assay(DESeq.ds)) 

### visualize rlog normalized data
#par(mfrow =c(2, 1)) # to plot images (in 2 rows, 1 col) 
# first, boxplots of non-transformed read counts
#boxplot (final.count.table, notch = True, main = "untransformed read counts", ylab = "read counts")

# next, box plots of rlog-transformed read counts
#boxplot (Final_rlog.counts, notch = TRUE, main = "rlog-transformed read counts", ylab = "rlog (read counts)")

### exploratory plots
###
# Pairwise correlation
?cor # cor () calculates the correlation between columns of a matrix
head(unfiltered_rlog.counts, n=10)
cor(unfiltered_rlog.counts[,1], unfiltered_rlog.counts[,5]) #default method is Pearson
cor(unfiltered_rlog.counts[,1], unfiltered_rlog.counts[,5], method = "spearman")

cor(unfiltered_rlog.counts, method = "spearman")


# clustering
# distance measure is d = 1 - r
distance.rlog <- as.dist (1 - cor(unfiltered_rlog.counts)) 
plot(hclust(distance.rlog)) 


#PCA
# Note that PCA (and clustering) should be done on normalized and preferably transformed read counts, so that the high variability of low read counts does not occlude informative data trends.
library (ggplot2)
plotPCA (unfiltered_DESeq.rlog) # Apply to DESeqDataSet; it may not be used with non-normalized read counts


########################################################################
# Sec 3: DGE analysis, result table and related graphs

# running the DGE analysis using the DESeq() function
# DESeq() function is basically a wrapper around the following three individual functions:
# (1) DESeq.ds <- estimateSizeFactors (DESeq.ds) # normalization between samples
# (2) DESeq.ds <- estimateDispersions (DESeq.ds) # gene-wise dispersion estimates 
# (3) DESeq.ds <- nbinomWaldTest (DESeq.ds) # this fits a negative binomial GLM (Logistic Regression) and applies Wald statistics to each gene
# The input object should contain raw read counts (non-normalized, untransformed), as the function will perform normalization and transformation under the hood; supplying anything but raw read counts will result in nonsensical results.

unfiltered.ds <-  DESeq(unfiltered.ds)
# additional info (#1-3 above) is added to the DESeqDataSet object; countData is not altered 

# generate a result table using the results() function, which extracts the base means across samples, log2 fold changes (log2FC), standard errors, test statistics etc.
unfilteredDE.results <- results (unfiltered.ds, contrast = c('condition','Tumor_sample','Non.tumor_sample'), alpha = 0.001) 
#contrast contains exactly three elements: the name of the "design" formula, the name of the numerator level for the fold change, and the name of the denominator level for the fold change
#alpha is the padj cutoff (padj in DESeq2 is FDR)
summary (unfilteredDE.results) 

#check the content of the DGE.results object (similar to a data frame)
head (unfilteredDE.results, 10)
table (unfilteredDE.results$padj < 0.001) #to build a contingency table of the counts at each factor level (here, TRUE/FALSE).
table (unfilteredDE.results$log2FoldChange > 2)
row.names (subset(unfilteredDE.results, padj < 0.001))

# shrink log2FoldChanges (for low expression) and generate an updated result table
unfilteredres.shrunk <- lfcShrink(unfiltered.ds, contrast = c('condition','Tumor_sample','Non.tumor_sample'), res=unfilteredDE.results, type = 'normal')
head(unfilteredres.shrunk, 10)# compare the result tables: head(res.shrunk) vs. head(DE.results)
class(unfilteredres.shrunk) #check the object type by class() function

### Now let's generate several exploratory plots of the DGE based on the result table
###
# MA plot
# y-axis: the expression change between conditions (log ratios, M); x-axis: the average expression strength of genes (average mean, 'A')
plotMA (unfilteredDE.results, alpha = 0.001, ylim = c(-5,5)) # genes that pass padj < 0.05 are colored in blue
plotMA (unfilteredres.shrunk, alpha = 0.001, ylim = c(-5,5))

# Volcano plot
# (1st time use only) Download and install the EnhancedVolcano package from Bioconductor
# BiocManager::install('EnhancedVolcano')

# load the package 
library(EnhancedVolcano)

EnhancedVolcano(unfilteredres.shrunk, lab = rownames(unfilteredres.shrunk), x ='log2FoldChange', y ='padj')
EnhancedVolcano(unfilteredres.shrunk, lab = rownames(unfilteredres.shrunk), x ='log2FoldChange', y ='padj', pCutoff = 10e-50, FCcutoff = 4) # the default pCutoff is 10e-6; default FCcutoff is abs(log2FC) > 1) 
# for more explanations and advanced features: https://www.bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
rab25 <- 'ENSG00000132698'
EnhancedVolcano(unfilteredres.shrunk, lab = rownames(unfilteredres.shrunk), x ='log2FoldChange', y ='padj', pCutoff = 10e-20, FCcutoff = 2, selectLab = rab25, drawConnectors = TRUE, labCol = 'black', labFace = 'bold', boxedLabels = TRUE)


# Histograms based on padj
# the frequency (y-axis) of certain values (x-axis) present in a data set
hist (unfilteredDE.results$padj)
hist (unfilteredDE.results$padj, 
      col='grey', breaks=50, 
      xlab = "adjusted p-value", ylab = "number of genes", main = "frequencies of adjusted p-values")


# Heatmap
## generating heatmaps using NMF::aheatmap()
#(first time use only) install the NMF package 
#BiocManager::install("NMF")

# load the library with the aheatmap() function
library (NMF) #load each time of use

# aheatmap() needs a matrix of values: in this example, a matrix of DE genes with the transformed read counts for each replicate
# identify genes with the desired adjusted p-value cut-off
unfilteredDEgenes <- row.names(subset(unfilteredDE.results, padj < 0.001))

# extract the normalized read counts (from Sec. 2) for DE genes into a matrix
#mat_DEgenes <- rlog.counts[DEgenes, ]
unfilteredmat_DEgenes <- unfiltered_rlog.counts[unfilteredDEgenes, c('Tumor_sample.1', 'Tumor_sample.2', 'Tumor_sample.3', 'Tumor_sample.4', 'Tumor_sample.5', 'Tumor_sample.6', 'Tumor_sample.7', 'Non.tumor_sample.1',  'Non.tumor_sample.2',  'Non.tumor_sample.3',  'Non.tumor_sample.4',  'Non.tumor_sample.5')]

# plot the normalized read counts of DE genes 
aheatmap (unfilteredmat_DEgenes, Rowv = NA , Colv = NA)

# combine the heatmap with hierarchical clustering (built-in of aheatmap)
dev.new()  # open a new plotting device with width and height of 10 inches
aheatmap(unfilteredmat_DEgenes,
         Rowv = TRUE, Colv = TRUE,  # add dendrograms to rows and columns
         distfun = "euclidean", hclustfun = "median")

# scale the read counts per gene to emphasize the sample-type-specific differences
dev.new()
aheatmap (unfilteredmat_DEgenes,
          Rowv = TRUE , Colv = TRUE ,
          distfun = "euclidean", hclustfun = "median",
          scale = "row") # values are transformed into distances from the center of the row-specific average: =(actual value - mean of the group) /(standard deviation); the colors now represent z-scores rather than the underlying read counts.

### generate output table/file
# export the results into a simple text file that can be opened with other programs for downstream analyses.
write.table (res.shrunk, file = "DESeq2results_D16.vs.D0.txt", sep = "\t", quote = FALSE)
write.table (subset(res.shrunk, padj<0.05), file = "DESeq2results_D16.vs.D0.significant2.txt", sep = "\t", quote = FALSE)

