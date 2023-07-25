setwd("C:/Users/albin/OneDrive/Desktop/SEM3/cancer genomics")

data<- read.table("GSE200146_counts.Cell.Lines.txt", 
                  header=TRUE, row.names=1, sep = '\t')
data= data[,-1]
head(data)


sizefactors= colSums(data)
columnsums = apply(data, 2, sum)


normalizedMatrix <- apply(data, 2, function(x)(x/sum(x))*1000000)
normalizedMatrixlog= log2(normalizedMatrix+1)
dim(normalizedMatrixlog)

# Assuming you have a numeric matrix named "data"
standardizedData <- matrix(NA, nrow = nrow(data), ncol = ncol(data))  # Create an empty matrix for standardized values
dim(standardizedData)

for (i in 1:100) {
  for (j in 1:100) {
    mean_value <- mean(normalizedMatrixlog[i, ])
    sd_value <- sd(normalizedMatrixlog[i, ])
    standardizedData[i, j] <- (normalizedMatrixlog[i, j] - mean_value) / sd_value
  }
}


mean_value <- mean(normalizedMatrixlog)
sd_value <- sd(normalizedMatrixlog)
z_scores <- (normalizedMatrixlog - mean_value) / sd_value

library(ComplexHeatmap)
library(circlize)
Heatmap(z_scores[1:100, ], col=colorRamp2(c(-2, 0, 2),c('orange', 'white', 'purple')))
saveRDS(z_scores, "z_scores.rds")

# Calculate variance for each gene
gene_var <- apply(data[1:100], 1, var)
print(gene_var)

#top 100genes w most variance
sorted_genes <- order(gene_var, decreasing = TRUE)
top_100_genes <- rownames(data)[sorted_genes[1:100]]
print(top_100_genes)

top100<-gene_var [1:100]
print(top100)
mean1_value <- mean(gene_var)
sd1_value <- sd(gene_var)
z1_scores <- (top100 - mean_value) / sd_value
View(z_scores)

#heatmap of zscore2
library(ComplexHeatmap)
library(circlize)
Heatmap(z_scores[1:100, ], col=colorRamp2(c(-2, 0, 2),c('orange', 'white', 'purple')))
