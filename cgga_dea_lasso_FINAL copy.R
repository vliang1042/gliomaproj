#install and load libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

if (!requireNamespace("BiocManager", force = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

install.packages("gprofiler2")

#loading necessary libraries 
library(limma)
library(tidyverse)
library(dplyr)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(factoextra)
library(glmnet)
library(pROC)
library(ggrepel)
library(ggplot2)
library(org.Hs.eg.db)
library(gprofiler2)
library(pheatmap)
library(caret)
library(data.table)
library(tableone)
library(survival)
library(survminer)

# Data downloadable from CGGA website. Change file path accordingly 

expression_file <- "~/Library/CloudStorage/OneDrive-UniversityofBirmingham/M07_Project/CGGA.mRNAseq_693.RSEM-genes.20200506.txt" 
cgga_expression_original <- read.table(expression_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

clinical_file <- "~/Library/CloudStorage/OneDrive-UniversityofBirmingham/M07_Project/CGGA.mRNAseq_693_clinical.20200506.txt"
cgga_clinical_original <- read.table(clinical_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Convert to numeric matrix
cgga_expression_original <- as.matrix(cgga_expression_original)
mode(cgga_expression_original) <- "numeric"

# keeping only WHO III and WHO IV and Primary tumour
cgga_clinical_original <- cgga_clinical_original %>%
  filter(
    Grade %in% c("WHO III", "WHO IV"),
    PRS_type == "Primary"
  )

# changing row names to CGGA ID 
rownames(cgga_clinical_original) <- cgga_clinical_original$CGGA_ID
cgga_clinical_original$CGGA_ID <- NULL 

# log2 transformation to reduce impact of extreme values 
cgga_expression_filtered <- log2(cgga_expression_original + 1)

# Filter low-expressed genes
threshold <- 1  # log2(FPKM + 1) expression threshold
min_samples <- round(0.2 * ncol(cgga_expression_filtered))  # must be expressed in at least 20% of samples

# Keep only genes expressed above threshold in enough samples
keep_genes <- rowSums(cgga_expression_filtered > threshold) >= min_samples
cgga_expression_filtered <- cgga_expression_filtered[keep_genes, ]

cat("Retained", sum(keep_genes), "genes out of", nrow(cgga_expression_original), "\n")

#keeping only grade variable from clinical data
cgga_grade <- cgga_clinical_original['Grade']
cgga_expression_filtered <- t(cgga_expression_filtered)

#filtering to match clinical data 
cgga_expression_filtered <- cgga_expression_filtered[rownames(cgga_expression_filtered) %in% rownames(cgga_grade), ]

#merged df for PCA 
merged_df <- merge(cgga_expression_filtered, cgga_grade, by = 'row.names')
rownames(merged_df) <- merged_df$Row.names
merged_df$Row.names <- NULL 

#checking for genes with 0 expression for all 
all_zero_cols <- sapply(cgga_expression_filtered, function(col) all(col == 0))
zero_expr_genes <- names(cgga_expression_filtered)[all_zero_cols]

write.csv(cgga_expression_filtered, "cgga_expression.csv")
write.csv(cgga_grade, "cgga_grade.csv")

cgga_expression_filtered <- t(cgga_expression_filtered)

#working out top 500 most variable genes 
gene_variances <- apply(cgga_expression_filtered, 1, var)
top_variable_genes <- names(sort(gene_variances, decreasing = TRUE)) [1:500]
expr_top <- cgga_expression_filtered[top_variable_genes, ]

#PCA 
expr_scaled <- t(scale(t(expr_top)))
pca_result <- prcomp(t(expr_scaled), scale. = FALSE)

group_factor <- factor(ifelse(cgga_clinical_original$Grade == "WHO III", 
                              "Lower Grade Glioma", 
                              "Glioblastoma"))

fviz_pca_ind(pca_result,
             geom.ind = "point",
             habillage = group_factor,  # This assigns color and creates legend
             addEllipses = TRUE,        # Optional: adds confidence ellipses
             title = "PCA of Top Variable Genes",
             legend.title = "Cancer Type")

# DEA

# load processed data 
cgga_clin <- read.csv("cgga_grade.csv")
cgga_expression <- read.csv("cgga_expression.csv")
meta_df <- read.delim("CGGA.mRNAseq_693_clinical.20200506.txt")

#make genes row names in clin and expression data
rownames(cgga_clin) <- cgga_clin$X
cgga_clin$X <- NULL 

rownames(cgga_expression) <- cgga_expression$X
cgga_expression$X <- NULL

cgga_expression <- t(cgga_expression)

# make grades factors 
cgga_clin$Grade <- factor(cgga_clin$Grade)
levels(cgga_clin$Grade) <- make.names(levels(cgga_clin$Grade))
cgga_clin$Grade <- relevel(cgga_clin$Grade, ref = "WHO.III")

# 2. Build the design matrix 
design <- model.matrix(~ Grade, data = cgga_clin)

# 3. Fit the linear model
fit <- lmFit(cgga_expression, design)

# 4. eBayes and extract results for WHO III and WHO IV 
fit2 <- eBayes(fit)

results <- topTable(fit2, adjust.method = "fdr", number = Inf)
head(results)

# Classify significance and direction
results <- results %>%
  mutate(
    regulation = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Up",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Get top 10 up and down regulated genes
top_genes <- results %>%
  filter(regulation != "NS") %>%
  arrange(desc(logFC)) %>%
  slice_head(n = 10) %>%
  bind_rows(
    results %>%
      filter(regulation != "NS") %>%
      arrange(logFC) %>%
      slice_head(n = 10)
  )

# Volcano plot with gene labels
ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = regulation)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "darkseagreen3", "NS" = "grey70"),
    breaks = c("Up", "Down"),
    name = "Regulation"
  ) +
  geom_text_repel(
    data = top_genes,
    aes(label = rownames(top_genes)), 
    color = "black",
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Differentially Expressed Genes in Glioblastoma v Lower Grade Glioma",
    x = "Log2 Fold Change (Glioblastoma v Lower Grade Glioma)",
    y = "-log10 Adjusted P-value (FDR)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# Enrichment Analysis and GSEA

## 1. Clean  IDs (if needed)
results_filtered <- results %>%
  mutate(SYMBOL = gsub("\\..*$", "", rownames(.)))

# Select only Up and Down regulated genes
deg_genes <- results_filtered %>%
  filter(regulation %in% c("Up", "Down")) %>%
  pull(SYMBOL)

## 2. GO enrichment
ego <- enrichGO(
  gene          = deg_genes,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff  = 0.05,
  readable      = TRUE,
  keyType       = "SYMBOL"
)

## 3. Barplot of top GO terms
barplot(ego, showCategory = 20, title = "GO Enrichment: Biological Processes") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),                  
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    plot.margin = ggplot2::margin(10, 10, 10, 20)          
  )

## 4. Prepare ranked gene list
gene_list <- results_filtered$logFC
names(gene_list) <- results_filtered$SYMBOL
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]

## 5. GSEA
gsea_result <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  keyType      = "SYMBOL",
  pvalueCutoff = 0.05
)

# Plot GSEA results (if available)
if (nrow(gsea_result@result) > 0) {
  ridgeplot(gsea_result, showCategory = 10) +
    ggtitle("Top Enriched GO Biological Processes") +
    xlab("Enrichment Score") + 
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 8),                    # reduce label size
      plot.title = element_text(hjust = 0.5, face = "bold"),   # center and bold the title
      plot.margin = ggplot2::margin(10, 10, 10, 10)            # spacing around the plot
    )
} else {
  message("No enriched pathways found. Consider relaxing filters or reviewing input genes.")
}

## 6. Display top enriched pathways as a table
if (nrow(gsea_result@result) > 0) {
  gsea_table <- gsea_result@result %>% 
    dplyr::select(ID, Description, pvalue, p.adjust, NES) %>% 
    dplyr::arrange(p.adjust) %>% 
    head(10)
  
  gsea_table$pvalue <- formatC(gsea_table$pvalue, format = "e", digits = 10)
  gsea_table$p.adjust <- formatC(gsea_table$p.adjust, format = "e", digits = 10)
  
  knitr::kable(gsea_table, caption = "Top Enriched GO Biological Processes")
} else {
  message("No pathways enriched at the specified cutoff.")
}

# overrepresentation analysis on up regulated genes 

# Filter for upregulated genes
up_genes <- results_filtered %>%
  filter(regulation == "Up")

## 2. GO enrichment
upreg_ego <- enrichGO(
  gene          = up_genes$SYMBOL,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff  = 0.05,
  readable      = TRUE,
  keyType       = "SYMBOL"
)

## 3. Barplot of top GO terms
barplot(upreg_ego, showCategory = 20, title = "GO Enrichment: Biological Processes in Upregulated Genes in GBM v LGG") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),                  
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    plot.margin = ggplot2::margin(10, 10, 10, 20)          
  )


# overrepresentation analysis on downregulated genes 

# Filter for down regulated genes
down_genes <- results_filtered %>%
  filter(regulation == "Down")

## 2. GO enrichment
downreg_ego <- enrichGO(
  gene          = down_genes$SYMBOL,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff  = 0.05,
  readable      = TRUE,
  keyType       = "SYMBOL"
)

## 3. Barplot of top GO terms
barplot(downreg_ego, showCategory = 20, title = "GO Enrichment: Biological Processes in Downregulated Genes in GBM v LGG") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),                  
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    plot.margin = ggplot2::margin(10, 10, 10, 20)          
  )

selected_genes_deg <- rownames(results_filtered[results_filtered$regulation %in% c("Up", "Down"), ])

labels <- cgga_clin$Grade

# split into train and test
set.seed(123)
train_idx <- sample(seq_len(ncol(cgga_expression)), size = 0.7 * ncol(cgga_expression))
train_expr <- cgga_expression[, train_idx]
test_expr  <- cgga_expression[, -train_idx]
train_labels <- labels[train_idx]
test_labels  <- labels[-train_idx]

# 3. Prepare feature matrices
train_data <- t(train_expr[selected_genes_deg, ])
test_data  <- t(test_expr[selected_genes_deg, ])

# Encode labels as binary: lower grade glioma = 0, glioblastoma = 1 
y_train <- ifelse(train_labels == "WHO.III", 0, 1)

# Fit logistic regression with LASSO using cross-validation
cv_fit <- cv.glmnet(
  x = as.matrix(train_data),
  y = y_train,
  family = "binomial",
  alpha = 1,  # LASSO
  nfolds = 5
)

# Fit final LASSO model at optimal lambda
lasso_model <- glmnet(train_data, y_train, alpha = 1, family = "binomial", lambda = cv_fit$lambda.min)

# Extract selected genes (non-zero coefficients)
coef_lasso <- coef(lasso_model)
selected_genes <- rownames(coef_lasso)[which(coef_lasso != 0)]
selected_genes <- selected_genes[!selected_genes %in% "(Intercept)"]

# Convert to a regular matrix and get names and values
coef_matrix <- as.matrix(coef_lasso)
nonzero_indices <- which(coef_matrix != 0)

# Get names and values of non-zero coefficients (excluding intercept)
nonzero_names <- rownames(coef_matrix)[nonzero_indices]
nonzero_values <- coef_matrix[nonzero_indices]

# Exclude intercept if present
keep <- nonzero_names != "(Intercept)"
selected_genes <- nonzero_names[keep]
selected_weights <- nonzero_values[keep]

# Combine names and weights into a data frame for easy viewing
lasso_results <- data.frame(Gene = selected_genes, Weight = selected_weights)
lasso_results <- lasso_results[order(-abs(lasso_results$Weight)), ]

cat("LASSO selected", length(selected_genes), "genes:\n")
print(selected_genes)

X_train_sel <- t(cgga_expression[selected_genes, train_idx])
X_test_sel  <- t(cgga_expression[selected_genes, -train_idx])

# Predict on test set
y_test <- ifelse(test_labels == "WHO.III", 0, 1)

# Fit logistic regression model
model_logit <- glm(y_train ~ ., data = data.frame(X_train_sel, y_train), family = binomial)

# Predict on test set
pred_probs <- predict(model_logit, newdata = data.frame(X_test_sel), type = "response")

# ROC Curve
roc_lr <- roc(response = y_test, predictor = as.numeric(pred_probs))

plot(roc_lr, print.auc = TRUE, col = "forestgreen", main = "Logistic Regression (LASSO) ROC Curve")

# Binomial Deviance vs Log(λ)
cv_lasso <- cv.glmnet(train_data, y_train, alpha = 1, family = "binomial")

plot(cv_lasso)
title("LASSO CV Plot (Binomial Deviance vs log(λ))", line = 2.5)

#Evaluation metrics 

# Convert test labels to binary (ensure they match predicted classes)
y_test_binary <- ifelse(y_test > 0.5, 1, 0)
pred_class <- ifelse(pred_probs > 0.5, 1, 0)


# Generate confusion matrix
library(caret)
conf_matrix <- confusionMatrix(as.factor(pred_class), as.factor(y_test_binary))
print(conf_matrix)

# Load ggplot2 for visualization
library(ggplot2)

# Extract confusion matrix table
cm_table <- as.data.frame(conf_matrix$table)

# Rename columns for clarity
colnames(cm_table) <- c("Predicted", "Actual", "Count")

# Plot confusion matrix as a heatmap
ggplot(data = cm_table, aes(x = Actual, y = Predicted, fill = Count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Count), vjust = 1.5, size = 6, color = "black") +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  labs(title = "Confusion Matrix", x = "Actual", y = "Predicted") +
  theme_minimal()

# Extract metrics
accuracy <- conf_matrix$overall["Accuracy"]
precision <- conf_matrix$byClass["Pos Pred Value"]     # Positive Predictive Value
recall <- conf_matrix$byClass["Sensitivity"]           # True Positive Rate
specificity <- conf_matrix$byClass["Specificity"]      # True Negative Rate

# Calculate F1-score
f1_score <- 2 * (precision * recall) / (precision + recall)

# Print the metrics
cat("Accuracy:", round(accuracy, 3), "\n")
cat("Precision:", round(precision, 3), "\n")
cat("Recall:", round(recall, 3), "\n")
cat("Specificity:", round(specificity, 3), "\n")
cat("F1 Score:", round(f1_score, 3), "\n")

# Visualisation 

# Volcano plot 
selected_genes_plot <- results[rownames(results) %in% selected_genes, ]


ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = regulation)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("Up" = "red", "Down" = "darkseagreen3", "NS" = "grey70"),
    breaks = c("Up", "Down"),
    name = "Regulation"
  ) +
  geom_text_repel(
    data = selected_genes_plot,
    aes(label = rownames(selected_genes_plot)),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    color = "black",
  ) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: LASSO-Selected Gene Labels",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted P-value"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# Heatmap of Selected Genes (LASSO)

# Subset metadata to only primary tumors
meta_df_primary <- meta_df[
  meta_df$PRS_type == "Primary" & meta_df$Grade %in% c("WHO III", "WHO IV"),
]

# Keep only samples with both expression and metadata
common_samples <- intersect(colnames(cgga_expression), meta_df_primary$CGGA_ID)

# Subset and reorder expression and metadata
expr_mat_final <- cgga_expression[, common_samples]
meta_df_final <- meta_df_primary[match(common_samples, meta_df_primary$CGGA_ID), ]

# Define condition
meta_df_final$condition <- ifelse(meta_df_final$Grade == "WHO IV", "GBM", "LGG")
meta_df_final$condition <- factor(meta_df_final$condition, levels = c("LGG", "GBM"))

# Subset expression matrix for selected genes and test samples
expr_heatmap <- expr_mat_final[selected_genes, -train_idx]

# Z-score scaling by gene (rows)
expr_scaled <- t(scale(t(expr_heatmap)))

# Prepare annotation for sample classes
annotation_col <- data.frame(Condition = y_test)
rownames(annotation_col) <- colnames(expr_scaled)

# Reorder columns (samples) by condition
sample_order <- order(y_test)  # y_test should be a factor with levels: LGG, GBM
expr_ordered <- expr_scaled[, sample_order]

# Update annotation to match new order

annotation_ordered <- annotation_col[sample_order, , drop = FALSE]

# Plot with fixed column order
pheatmap(
  expr_ordered,
  annotation_col = annotation_ordered,
  cluster_rows = TRUE,       # Keep clustering for genes
  cluster_cols = FALSE,      # Disable clustering for samples
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  border_color = "grey30",
  main = "Heatmap of LASSO Genes\nSamples Ordered by LGG vs GBM"
)

boxplot(as.numeric(expr_heatmap["MMP9", ]) ~ y_test,
        main = "Expression of MMP9",
        ylab = "log2(FPKM + 1) Expression", xlab = "Condition", 
        names = c("Lower Grade Glioma", "Glioblastoma"), 
        col = c("blue", "purple"))

boxplot(as.numeric(expr_heatmap["ADM", ]) ~ y_test,
        main = "Expression of ADM",
        ylab = "log2(FPKM + 1) Expression", xlab = "Condition", 
        names = c("Lower Grade Glioma", "Glioblastoma"),
        col = c("blue", "purple"))

boxplot(as.numeric(expr_heatmap["KLRK1", ]) ~ y_test,
        main = "Expression of KLRK1",
        ylab = "log2(FPKM + 1) Expression", xlab = "Condition", 
        names = c("Lower Grade Glioma", "Glioblastoma"), 
        col = c("blue", "purple"))

boxplot(as.numeric(expr_heatmap["NOX4", ]) ~ y_test,
        main = "Expression of NOX4",
        ylab = "log2(FPKM + 1) Expression", xlab = "Condition", 
        names = c("Lower Grade Glioma", "Glioblastoma"), 
        col = c("blue", "purple"))

boxplot(as.numeric(expr_heatmap["IL13RA2", ]) ~ y_test,
        main = "Expression of IL13RA2",
        ylab = "log2(FPKM + 1) Expression", xlab = "Condition", 
        names = c("Lower Grade Glioma", "Glioblastoma"), 
        col = c("blue", "purple"))

library(ggplot2)

# Prepare data
df <- data.frame(
  expression = as.numeric(expr_heatmap["IL13RA2", ]),
  condition = factor(y_test, labels = c("Lower Grade Glioma", "Glioblastoma"))
)

# Create ggplot boxplot
ggplot(df, aes(x = condition, y = expression, fill = condition)) +
  geom_boxplot() +
  labs(
    title = "Expression of IL13RA2",
    x = "Condition",
    y = "log2(FPKM + 1) Expression"
  ) +
  theme_minimal()

