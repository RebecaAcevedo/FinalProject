### Generate our HPC raw count (filtered) table with unique row names

finalcounts <- read.table('counts_us.txt', sep="\t", header=TRUE)
finalcounts$Geneid <- gsub("\\..*", "", finalcounts$Geneid)

finalcounts2 <- subset(finalcounts, select=c(7:18))
View(finalcounts2)

final.count.table <- data.frame(finalcounts2)
View(final.count.table)
names(final.count.table)

colnames(final.count.table)[colnames(final.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.SRR349741_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.1"
colnames(final.count.table)[colnames(final.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.SRR349742_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.2"
colnames(final.count.table)[colnames(final.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.SRR349743_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.3"
colnames(final.count.table)[colnames(final.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.SRR349744_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.4"
colnames(final.count.table)[colnames(final.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.SRR349745_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.5"
colnames(final.count.table)[colnames(final.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.SRR349746_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.6"
colnames(final.count.table)[colnames(final.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.SRR349747_Aligned.sortedByCoord.out.bam"] ="Tumor_sample.7"
colnames(final.count.table)[colnames(final.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.SRR349748_Aligned.sortedByCoord.out.bam"] ="Non.tumor_sample.1"
colnames(final.count.table)[colnames(final.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.SRR349749_Aligned.sortedByCoord.out.bam"] ="Non.tumor_sample.2"
colnames(final.count.table)[colnames(final.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.SRR349750_Aligned.sortedByCoord.out.bam"] ="Non.tumor_sample.3"
colnames(final.count.table)[colnames(final.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.SRR349751_Aligned.sortedByCoord.out.bam"] ="Non.tumor_sample.4"
colnames(final.count.table)[colnames(final.count.table) == "X.xdisk.guangyao.416a516a.Group15.chance.STAR.SRR349752_Aligned.sortedByCoord.out.bam"] ="Non.tumor_sample.5"

gene_id_unique <- finalcounts$Geneid
gene_id_unique <- make.names(gene_id_unique, unique=TRUE)
row.names(final.count.table) <- gene_id_unique
View(final.count.table)
dim(final.count.table)


### Generate recount2 raw count table with unique row names

recountgenes <- read.table("counts_gene.tsv", header=TRUE)
recountgenes$gene_id <- gsub("\\..*", "", recountgenes$gene_id)
#write.table(recountgenes, file = "output_them.txt", sep = "\t", quote = FALSE)

recounts2 <- subset(recountgenes, select=c(1:12))
recount.count.table <- data.frame(recounts2)

colnames(recount.count.table)[1:12] <- c("Tumor_sample.1", "Tumor_sample.2", "Tumor_sample.3", "Tumor_sample.4", "Tumor_sample.5", "Tumor_sample.6", "Tumor_sample.7", "Non.tumor_sample.1", "Non.tumor_sample.2", "Non.tumor_sample.3", "Non.tumor_sample.4", "Non.tumor_sample.5")
names(recount.count.table)

gene_id_unique <- recountgenes$gene_id
gene_id_unique <- make.names(gene_id_unique, unique=TRUE)
row.names(recount.count.table) <- gene_id_unique
View(recount.count.table)
dim(recount.count.table)


### merge count tables (first 12 columns our counts filtered/ second 12 are recount2 counts)

  # raw counts
merged.tables <- merge(final.count.table, recount.count.table, by = "row.names", no.dups = TRUE)
View(merged.tables)
dim(merged.tables)

rownames(merged.tables) <- merged.tables[,1]
merged.tables <- merged.tables[,-1]
dim(merged.tables)

  # rlog normalized counts
merged.rlog.counts <- merge(Final_rlog.counts, RecountFinal_rlog.counts, by = "row.names", no.dups = TRUE)
View(merged.rlog.counts)
dim(merged.rlog.counts)

rownames(merged.rlog.counts) <- merged.rlog.counts[,1]
merged.rlog.counts <- merged.rlog.counts[,-1]
dim(merged.tables)

### correlation analysis
  #raw counts 

library(ggplot2)

#SRR349741
dev.new()
par(mar=c(7, 5, 4, 2))
cor1 <- cor(merged.tables[,1], merged.tables[,13])
cor1 <- round(cor1, digits = 4)
plot(merged.tables[,1], merged.tables[,13], main = paste('Correlation = ', cor1), xlab = "Tumor1 HPC raw counts", ylab = "Tumor1 Recount raw counts")

#SRR349742
dev.new()
par(mar=c(7, 5, 4, 2))
cor2 <- cor(merged.tables[,2], merged.tables[,14])
cor2 <- round(cor2, digits = 4)
plot(merged.tables[,2], merged.tables[,14], main = paste('Correlation = ', cor2), xlab = "Tumor2 HPC raw counts", ylab = "Tumor2 Recount raw counts")

#SRR349743
dev.new()
par(mar=c(7, 5, 4, 2))
cor3 <- cor(merged.tables[,3], merged.tables[,15])
cor3 <- round(cor3, digits = 4)
plot(merged.tables[,3], merged.tables[,15], main = paste('Correlation = ', cor3), xlab = "Tumor3 HPC raw counts", ylab = "Tumor3 Recount raw counts")

#SRR349744
dev.new()
par(mar=c(7, 5, 4, 2))
cor4 <- cor(merged.tables[,4], merged.tables[,16])
cor4 <- round(cor4, digits = 4)
plot(merged.tables[,4], merged.tables[,16], main = paste('Correlation = ', cor4), xlab = "Tumor4 HPC raw counts", ylab = "Tumor4 Recount raw counts")

#SRR349745
dev.new()
par(mar=c(7, 5, 4, 2))
cor5 <- cor(merged.tables[,5], merged.tables[,17])
cor5 <- round(cor5, digits = 4)
plot(merged.tables[,5], merged.tables[,17], main = paste('Correlation = ', cor5), xlab = "Tumor5 HPC raw counts", ylab = "Tumor5 Recount raw counts")

#SRR349746
dev.new()
par(mar=c(7, 5, 4, 2))
cor6 <- cor(merged.tables[,6], merged.tables[,18])
cor6 <- round(cor6, digits = 4)
plot(merged.tables[,6], merged.tables[,18], main = paste('Correlation = ', cor6), xlab = "Tumor6 HPC raw counts", ylab = "Tumor6 Recount raw counts")

#SRR349747
dev.new()
par(mar=c(7, 5, 4, 2))
cor7 <- cor(merged.tables[,7], merged.tables[,19])
cor7 <- round(cor7, digits = 4)
plot(merged.tables[,7], merged.tables[,19], main = paste('Correlation = ', cor7), xlab = "Tumor7 HPC raw counts", ylab = "Tumor7 Recount raw counts")

#SRR349748
dev.new()
par(mar=c(7, 5, 4, 2))
cor8 <- cor(merged.tables[,8], merged.tables[,20])
cor8 <- round(cor8, digits = 4)
plot(merged.tables[,8], merged.tables[,20], main = paste('Correlation = ', cor8), xlab = "NonTumor1 HPC raw counts", ylab = "NonTumor1 Recount raw counts")

#SRR349749
dev.new()
par(mar=c(7, 5, 4, 2))
cor9 <- cor(merged.tables[,9], merged.tables[,21])
cor9 <- round(cor9, digits = 4)
plot(merged.tables[,9], merged.tables[,21], main = paste('Correlation = ', cor9), xlab = "NonTumor2 HPC raw counts", ylab = "NonTumor2 Recount raw counts")

#SRR349750
dev.new()
par(mar=c(7, 5, 4, 2))
cor10 <- cor(merged.tables[,10], merged.tables[,22])
cor10 <- round(cor10, digits = 4)
plot(merged.tables[,10], merged.tables[,22], main = paste('Correlation = ', cor10), xlab = "NonTumor3 HPC raw counts", ylab = "NonTumor3 Recount raw counts")

#SRR349751
dev.new()
par(mar=c(7, 5, 4, 2))
cor11 <- cor(merged.tables[,11], merged.tables[,23])
cor11 <- round(cor11, digits = 4)
plot(merged.tables[,11], merged.tables[,23], main = paste('Correlation = ', cor11), xlab = "NonTumor4 HPC raw counts", ylab = "NonTumor4 Recount raw counts")

#SRR349752
dev.new()
par(mar=c(7, 5, 4, 2))
cor12 <- cor(merged.tables[,12], merged.tables[,24])
cor12 <- round(cor12, digits = 4)
plot(merged.tables[,12], merged.tables[,24], main = paste('Correlation = ', cor12), xlab = "NonTumor5 HPC raw counts", ylab = "NonTumor5 Recount raw counts")



### unfiltered ######################################################################

### Generate our HPC raw count(unfiltered) table with unique row names

unfilteredcounts <- read.table('counts_unfiltered.txt', sep="\t", header=TRUE)
unfilteredcounts$Geneid <- gsub("\\..*", "", unfilteredcounts$Geneid)

unfilteredcounts2 <- subset(unfilteredcounts, select=c(7:18))
View(unfilteredcounts2)

unfiltered.count.table <- data.frame(unfilteredcounts2)
View(unfiltered.count.table)
names(unfiltered.count.table)

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

gene_id_unique <- unfilteredcounts$Geneid
gene_id_unique <- make.names(gene_id_unique, unique=TRUE)
row.names(unfiltered.count.table) <- gene_id_unique
View(unfiltered.count.table)
dim(unfiltered.count.table)


### Generate recount2 raw count table with unique row names

recountgenes <- read.table("counts_gene.tsv", header=TRUE)
recountgenes$gene_id <- gsub("\\..*", "", recountgenes$gene_id)
#write.table(recountgenes, file = "output_them.txt", sep = "\t", quote = FALSE)

recounts2 <- subset(recountgenes, select=c(1:12))
recount.count.table <- data.frame(recounts2)

colnames(recount.count.table)[1:12] <- c("Tumor_sample.1", "Tumor_sample.2", "Tumor_sample.3", "Tumor_sample.4", "Tumor_sample.5", "Tumor_sample.6", "Tumor_sample.7", "Non.tumor_sample.1", "Non.tumor_sample.2", "Non.tumor_sample.3", "Non.tumor_sample.4", "Non.tumor_sample.5")
names(recount.count.table)

gene_id_unique <- recountgenes$gene_id
gene_id_unique <- make.names(gene_id_unique, unique=TRUE)
row.names(recount.count.table) <- gene_id_unique
View(recount.count.table)
dim(recount.count.table)


### merge count tables (first 12 columns our counts/ second 12 are recount2 counts)

# raw counts unfiltered
merged.tables <- merge(unfiltered.count.table, recount.count.table, by = "row.names", no.dups = TRUE)
#write.table (merged.tables, file = "raw unfiltered_recount merged.txt", sep = "\t", quote = FALSE)
View(merged.tables)
dim(merged.tables)

rownames(merged.tables) <- merged.tables[,1]
merged.tables <- merged.tables[,-1]
dim(merged.tables)



# plot correlation

library(ggplot2)

#SRR349741
dev.new()
par(mar=c(7, 5, 4, 2))
cor1 <- cor(merged.tables[,1], merged.tables[,13])
cor1 <- round(cor1, digits = 4)
plot(merged.tables[,1], merged.tables[,13], main = paste('Correlation = ', cor1), xlab = "Tumor1 HPC raw counts (unfiltered)", ylab = "Tumor1 Recount raw counts")


#SRR349742
dev.new()
par(mar=c(7, 5, 4, 2))
cor2 <- cor(merged.tables[,2], merged.tables[,14])
cor2 <- round(cor2, digits = 4)
plot(merged.tables[,2], merged.tables[,14], main = paste('Correlation = ', cor2), xlab = "Tumor2 HPC raw counts (unfiltered)", ylab = "Tumor2 Recount raw counts")

#SRR349743
dev.new()
par(mar=c(7, 5, 4, 2))
cor3 <- cor(merged.tables[,3], merged.tables[,15])
cor3 <- round(cor3, digits = 4)
plot(merged.tables[,3], merged.tables[,15], main = paste('Correlation = ', cor3), xlab = "Tumor3 HPC raw counts (unfiltered)", ylab = "Tumor3 Recount raw counts")

#SRR349744
dev.new()
par(mar=c(7, 5, 4, 2))
cor4 <- cor(merged.tables[,4], merged.tables[,16])
cor4 <- round(cor4, digits = 4)
plot(merged.tables[,4], merged.tables[,16], main = paste('Correlation = ', cor4), xlab = "Tumor4 HPC raw counts (unfiltered)", ylab = "Tumor4 Recount raw counts")

#SRR349745
dev.new()
par(mar=c(7, 5, 4, 2))
cor5 <- cor(merged.tables[,5], merged.tables[,17])
cor5 <- round(cor5, digits = 4)
plot(merged.tables[,5], merged.tables[,17], main = paste('Correlation = ', cor5), xlab = "Tumor5 HPC raw counts (unfiltered)", ylab = "Tumor5 Recount raw counts")

#SRR349746
dev.new()
par(mar=c(7, 5, 4, 2))
cor6 <- cor(merged.tables[,6], merged.tables[,18])
cor6 <- round(cor6, digits = 4)
plot(merged.tables[,6], merged.tables[,18], main = paste('Correlation = ', cor6), xlab = "Tumor6 HPC raw counts (unfiltered)", ylab = "Tumor6 Recount raw counts")

#SRR349747
dev.new()
par(mar=c(7, 5, 4, 2))
cor7 <- cor(merged.tables[,7], merged.tables[,19])
cor7 <- round(cor7, digits = 4)
plot(merged.tables[,7], merged.tables[,19], main = paste('Correlation = ', cor7), xlab = "Tumor7 HPC raw counts (unfiltered)", ylab = "Tumor7 Recount raw counts")

#SRR349748
dev.new()
par(mar=c(7, 5, 4, 2))
cor8 <- cor(merged.tables[,8], merged.tables[,20])
cor8 <- round(cor8, digits = 4)
plot(merged.tables[,8], merged.tables[,20], main = paste('Correlation = ', cor8), xlab = "NonTumor1 HPC raw counts (unfiltered)", ylab = "NonTumor1 Recount raw counts")

#SRR349749
dev.new()
par(mar=c(7, 5, 4, 2))
cor9 <- cor(merged.tables[,9], merged.tables[,21])
cor9 <- round(cor9, digits = 4)
plot(merged.tables[,9], merged.tables[,21], main = paste('Correlation = ', cor9), xlab = "NonTumor2 HPC raw counts (unfiltered)", ylab = "NonTumor2 Recount raw counts")

#SRR349750
dev.new()
par(mar=c(7, 5, 4, 2))
cor10 <- cor(merged.tables[,10], merged.tables[,22])
cor10 <- round(cor10, digits = 4)
plot(merged.tables[,10], merged.tables[,22], main = paste('Correlation = ', cor10), xlab = "NonTumor3 HPC raw counts (unfiltered)", ylab = "NonTumor3 Recount raw counts")

#SRR349751
dev.new()
par(mar=c(7, 5, 4, 2))
cor11 <- cor(merged.tables[,11], merged.tables[,23])
cor11 <- round(cor11, digits = 4)
plot(merged.tables[,11], merged.tables[,23], main = paste('Correlation = ', cor11), xlab = "NonTumor4 HPC raw counts (unfiltered)", ylab = "NonTumor4 Recount raw counts")

#SRR349752
dev.new()
par(mar=c(7, 5, 4, 2))
cor12 <- cor(merged.tables[,12], merged.tables[,24])
cor12 <- round(cor12, digits = 4)
plot(merged.tables[,12], merged.tables[,24], main = paste('Correlation = ', cor12), xlab = "NonTumor5 HPC raw counts (unfiltered)", ylab = "NonTumor5 Recount raw counts")