library(DESeq2)

input_file <- commandArgs(trailingOnly = TRUE)[1]
output_file <- commandArgs(trailingOnly = TRUE)[2]

print(input_file)
countdata <- read.table(input_file, header=TRUE, row.names=1, fill = TRUE)
head(countdata)
countdata <- countdata[, 6:ncol(countdata)] # Remove unwanted columns
countdata <- as.matrix(countdata)
samples <- colnames(countdata)
(coldata <- data.frame(row.names=samples, sampleName=samples))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~1)
dds <- DESeq(dds)
SizeFactorsInput <- 1 / sizeFactors(dds)

normalized_counts <- counts(dds, normalized=TRUE)
normframe <- data.frame(normalized_counts)
normoverone <- normframe/do.call(pmax, normframe)
write.table(normalized_counts, file="results/NormalizedCounts.csv", sep=",", quote=F, col.names=NA)
write.table(normoverone, file="results/NormOverOne.csv", sep=",", quote=F, col.names=NA)
write.table(SizeFactorsInput, file=output_file, sep="\t", quote=FALSE, col.names=FALSE)


