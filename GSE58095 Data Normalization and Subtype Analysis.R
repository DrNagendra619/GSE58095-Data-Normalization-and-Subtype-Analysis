# 1. Load Libraries for GSE58095 Analysis
library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)

# 2. Load Dataset (GSE58095 - Systemic Sclerosis)
# This uses the Illumina HumanHT-12 V4.0 platform.
gse <- getGEO("GSE58095", GSEMatrix = TRUE)
metadata <- pData(gse[[1]])
expression_data <- exprs(gse[[1]])

# 3. Check for Missing Values in GSE58095
cat("Missing values in GSE58095 expression data:", sum(is.na(expression_data)), "\n")

# 4. Boxplot Before Normalization (GSE58095)
png(file = file.path("GSE58095_boxplot_before_normalization.png"))
boxplot(expression_data, main = "GSE58095 (Illumina) - Before Normalization", las = 2, outline = FALSE)
dev.off()

# 5. Apply Normalization to GSE58095
expression_data <- normalizeBetweenArrays(expression_data, method = "quantile")

# 6. Boxplot After Normalization (GSE58095)
png(file = file.path("GSE58095_boxplot_after_normalization.png"))
boxplot(expression_data, main = "GSE58095 (Illumina) - After Normalization", las = 2, outline = FALSE)
dev.off()

# 7. Filter Low-Expressed Genes from GSE58095
gene_means <- rowMeans(expression_data)
threshold <- quantile(gene_means, probs = 0.25)
expression_data <- expression_data[gene_means > threshold, ]
cat("GSE58095 Dimensions after filtering:", dim(expression_data), "\n")

# 8. Extract Subtypes (Healthy vs Diffuse vs Limited)
# We combine all metadata columns into text to ensure we find the disease terms
# even though the titles are just 'PATIENT01'.
meta_text <- apply(metadata, 1, paste, collapse = " ")

subtype <- dplyr::case_when(
  grepl("Control", meta_text, ignore.case = TRUE) ~ "Healthy Control",
  grepl("diffuse", meta_text, ignore.case = TRUE) ~ "Diffuse Scleroderma",
  grepl("limited", meta_text, ignore.case = TRUE) ~ "Limited Scleroderma",
  TRUE ~ "Unknown"
)
subtype <- as.factor(subtype)

cat("GSE58095 Subtype distribution:\n")
table(subtype)

# 9. Generate PCA Plot for GSE58095
pca_result <- prcomp(t(expression_data), scale. = TRUE)
pca_df <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], Subtype = subtype)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subtype)) +
  geom_point(size = 3, alpha = 0.8) +
  ggtitle("GSE58095 PCA: Scleroderma Subtypes (Illumina Array)") +
  theme_minimal()

ggsave(file.path("GSE58095_pca_plot.png"), plot = pca_plot)

# 10. Save Processed Data for GSE58095
save(expression_data, metadata, subtype, file = file.path("GSE58095_processed_data_step1.RData"))
