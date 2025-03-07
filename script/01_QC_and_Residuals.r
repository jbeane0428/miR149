
# Setup working environment
library(Biobase)
library(SummarizedExperiment)
library(limma)
library(edgeR)

# ==============================================================================
# Calculate residual expression
# ==============================================================================
count_bx <- readRDS("../input/PCGA_miRNA_raw_count_table_BX.rds")
dim(count_bx)

index <- readRDS("../input/PCGA_miRNA_Annotation_GEO_upload.rds")
index_bx <- index[index$Sample_Type == "Biopsy" & index$Filter == "Included", ]

count_bx <- count_bx[, row.names(index_bx)]

# ==============================================================================
# Calculate residual
# ==============================================================================
data.dge <- resid <- c()
data.dge <- DGEList(counts=count_bx, lib.size = as.numeric(index_bx[colnames(count_bx),]$Total_Aligned_by_mirdeep2))
logcpm_bx <- voom(data.dge, normalize.method="quantile")$E
keep <- row.names(logcpm_bx[rowSums(logcpm_bx >= 1 ) > nrow(index_bx)/2, ])
data.dge <- data.dge[keep, keep.lib.size=TRUE]
design.r <- model.matrix(~as.factor(index_bx$Batch))
v.r <- voom(data.dge,design=design.r,normalize.method="quantile")
fit.r <- lmFit(v.r,design.r)
fit.r <- eBayes(fit.r)
resid <- residuals(fit.r, v.r)

saveRDS(resid, "../output/PCGA_miRNA_Residual_BX.rds")
saveRDS(index_bx, "../output/PCGA_miRNA_index_BX.rds")

# ==============================================================================
# Filtered down gene expression and 
# ==============================================================================
# Get gene residuals
mrna_residual <- readRDS("../input/resid_bx.rds")
mrna_residual <- mrna_residual[,colnames(mrna_residual) %in% index_bx$RNA_ID]
mrna_residual <- mrna_residual[, as.character(index_bx$RNA_ID)]
colnames(mrna_residual) <- row.names(index_bx)
saveRDS(mrna_residual, "../output/PCGA_mrna_residual_BX.rds")

module_scores <- index_bx[,grepl("Module_", names(index_bx))]
names(module_scores) <- gsub("_", "", names(module_scores))
saveRDS(module_scores, "../output/PCGA_mRNA_modules_score_BX.rds")

