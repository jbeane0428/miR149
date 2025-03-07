
# Setup working environment
library(Biobase)
library(tidyverse)
library(SummarizedExperiment)
library(ggplot2)
library(gplots)
library(limma)
library(edgeR)
library(wesanderson)
library(ggsignif)
library(reshape2)
library(corrplot)
library(RColorBrewer)

# ==============================================================================
# Necessary functions
# ==========================================================================
cor2pvalue = function(r, n) {
	t <- (r*sqrt(n-2))/sqrt(1-r^2)
	p <- 2*(1 - pt(abs(t),(n-2)))
	return(p)
}

# ==============================================================================
# DE Analysis between progressive and regressive of Proliferative subtype
# ==============================================================================
#
# Module 9 related miRNAs
# ------------------------------------------------------------------------------
gene_module <- readRDS("../input/PCGA_Gene_Module.rds")
mir.module.gs <- readRDS("../output/miRNA_by_PCGA_Gene_Module.rds")

mirna <- readRDS("../output/PCGA_miRNA_Residual_BX.rds")
index_bx <- readRDS("../output/PCGA_miRNA_index_BX.rds")

# Select progressive vs. regressive within Proliferative subtype
index_bx <- index_bx[index_bx$Molecular_Subtype == "Proliferative", ]
index_bx <- index_bx[index_bx$Progression_Status %in% c("Progressive/Persistent", "Regressive"),]
index_bx[["Progression_Status"]] <- factor(index_bx[["Progression_Status"]], levels = c("Progressive/Persistent", "Regressive"))
mirna_limma <- mirna_res[, colnames(mirna_res) %in% row.names(index_bx)]

design.bx <- model.matrix(~index_bx$Progression_Status)
pat.block.bx <- as.factor(as.character(index_bx$Patient))
dupcor.bx <- duplicateCorrelation(mirna_limma, design.bx, block=pat.block.bx)
dupcor.bx$consensus.correlation 
fit.bx <- lmFit(mirna_limma, design.bx, block=pat.block.bx, correlation=dupcor.bx$consensus.correlation)
fit.bx <- eBayes(fit.bx)
res_bx<-topTable(fit.bx,coef=2, n=nrow(mirna_limma))
res_bx <- as.data.frame(res_bx)
res_bx[mir.module.gs[["Module9"]], ]

#
# Target genes of miR-149-5p
# ------------------------------------------------------------------------------
fil.net <- readRDS("../output/PCGA_BX_Filtered_miRNA_Gene_Network.rds")
plot.net <- fil.net["hsa-miR-149-5p", colnames(fil.net) %in% row.names(gene_module[gene_module$Gene_Module_Number == "Module9", ])]
mir.target <- names(plot.net)[plot.net == 1]

nlrc5.target <- c("ENSG00000206503", "ENSG00000234745", "ENSG00000204525",
            "ENSG00000166710", "ENSG00000204264", "ENSG00000240065", "ENSG00000168394")

mrna <- readRDS("../input/resid_bx.rds")
mrna_index <- readRDS("../input/samples_annotation_v3.rds")
mrna_index <- mrna_index[row.names(mrna_index) %in% colnames(mrna), c("Patient", "Molecular_Subtype", "Progression_Status")]

index_mrna_bx <- mrna_index[mrna_index$Molecular_Subtype=="Proliferative", ]
index_mrna_bx <- index_mrna_bx[index_mrna_bx$Progression_Status %in% c("Progressive/Persistent", "Regressive"),]
index_mrna_bx[["Progression_Status"]] <- factor(index_mrna_bx[["Progression_Status"]], levels = c("Progressive/Persistent", "Regressive"))
mrna_limma <- mrna[, row.names(index_mrna_bx)]

design.bx <- model.matrix(~index_mrna_bx$Progression_Status)
pat.block.bx <- as.factor(as.character(index_mrna_bx$Patient))
dupcor.bx <- duplicateCorrelation(mrna_limma, design.bx, block=pat.block.bx)
dupcor.bx$consensus.correlation 
fit.bx <- lmFit(mrna_limma, design.bx, block=pat.block.bx, correlation=dupcor.bx$consensus.correlation)
fit.bx <- eBayes(fit.bx)
res_bx<-topTable(fit.bx,coef=2, n=nrow(mrna_limma))
res_bx <- as.data.frame(res_bx)

res_bx[mir.target,]

res_bx[nlrc5.target,]

# ==============================================================================
# Plots
# ==============================================================================
progress_pal <- wes_palette("Chevalier1")[c(1,3,2,4)]

#
# Figure 2D
# ------------------------------------------------------------------------------
mir.mat <- mirna_limma[mir.module.gs[["Module9"]], ]
out.dat <- cbind(index_bx, t(mir.mat))

# Select progressive vs. regressive within Proliferative subtype
out.dat <- out.dat[out.dat$Molecular_Subtype == "Proliferative", ]
out.dat[["Progression_Status"]] <- factor(out.dat[["Progression_Status"]], levels = c("Progressive/Persistent", "Regressive", "Normal/Stable", "UNK"))
out.dat <- out.dat[out.dat$Progression_Status %in% c("Progressive/Persistent", "Regressive"),]

plot.mat <- out.dat[,c("Progression_Status", mir.module.gs[["Module9"]])]
out.table <- gather(plot.mat, Name, Value, colnames(plot.mat)[2]:colnames(plot.mat)[ncol(plot.mat)])
levels(out.table$Progression_Status) <- c("Progressive/\n  persistent\n     N=14",
										"Regresssive\n     N=14", NA, NA)
out.table$Progression_Status <- as.character(out.table$Progression_Status)

pdf("../output/Figure_2/Figure_2D.pdf", width=4, height=4)
    ggplot(out.table, aes(x=Name, y=Value, fill=Progression_Status)) + 
		geom_point(aes(color=Progression_Status), shape=15, size=0, alpha=0) +
		guides(color = guide_legend(override.aes = list(size = 5, alpha=1))) +
        geom_boxplot(outlier.alpha = 0) + 
		guides(fill = FALSE) +
		geom_point(position=position_jitterdodge(jitter.width=0.2), size=0.6) + 
		geom_signif(
			y_position = c(2.4), xmin = c(3.8), xmax = c(4.2),
			annotation = c("*"), tip_length = 0.01
		) +
		ylim(-2.5, 2.5) +
        scale_fill_manual(values=progress_pal) +
        scale_color_manual(values=progress_pal) +
        labs(x="Module 9 miRNA", y = "Residual Exprsesion") +
        theme_classic() +
        theme(axis.text.x=element_text(angle = -45, hjust = 0, size = 10),
            legend.position="bottom", 
			legend.title = element_blank(),
			panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank())
dev.off()

#
# Figure 3A
# ------------------------------------------------------------------------------
# Format data
mrna.mat <- mrna_limma[mir.target, ]
row.names(mrna.mat) <- gene_module[row.names(mrna.mat), ]$hgnc_symbol
out.dat <- cbind(index_mrna_bx, t(mrna.mat))

plot.mat <- out.dat[,c("Progression_Status", row.names(mrna.mat))]
out.table <- gather(plot.mat, Name, Value, colnames(plot.mat)[2]:colnames(plot.mat)[ncol(plot.mat)])
levels(out.table$Progression_Status) <- c("Progressive/\n  persistent\n     N=15",
										"Regresssive\n     N=15", NA, NA)
out.table$Progression_Status <- as.character(out.table$Progression_Status)

pdf("../output/Figure_3/Figure_3A.pdf", width=8, height=6)
    ggplot(out.table, aes(x=Name, y=Value, fill=Progression_Status)) + 
		geom_point(aes(color=Progression_Status), shape=15, size=0, alpha=0) +
		guides(color = guide_legend(override.aes = list(size = 5, alpha=1))) +
        geom_boxplot(outlier.alpha = 0) + 
		guides(fill = FALSE) +
		geom_point(position=position_jitterdodge(jitter.width=0.2)) + 
		geom_signif(
			y_position = c(2.6, 2.6, 2.6), xmin = c(2.7, 3.7, 4.7), xmax = c(3.3, 4.3, 5.3),
			annotation = c("*", "*", "*"), tip_length = 0.01
		) +
		ylim(-2, 3) +
        scale_fill_manual(values=progress_pal) +
        scale_color_manual(values=progress_pal) +
        labs(x="miR-149-5p Targets in Module 9", y = "Residual Exprsesion") +
        theme_classic() +
        theme(axis.text.x=element_text(angle = -45, hjust = 0, size = 10),
            legend.position="bottom", 
			legend.title = element_blank(),
			panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank())
dev.off()

#
# Figure 3B
# ------------------------------------------------------------------------------
mrna.mat <- mrna_limma[nlrc5.target, ]
row.names(mrna.mat) <- gene_module[row.names(mrna.mat), ]$hgnc_symbol
out.dat <- cbind(index_mrna_bx, t(mrna.mat))

plot.mat <- out.dat[,c("Progression_Status", row.names(mrna.mat))]
out.table <- gather(plot.mat, Name, Value, colnames(plot.mat)[2]:colnames(plot.mat)[ncol(plot.mat)])
levels(out.table$Progression_Status) <- c("Progressive/\n  persistent\n     N=15",
										"Regresssive\n     N=15", NA, NA)
out.table$Progression_Status <- as.character(out.table$Progression_Status)

pdf("../output/Figure_3/Figure_3B.pdf", width=8, height=6)
    ggplot(out.table, aes(x=Name, y=Value, fill=Progression_Status)) + 
		geom_point(aes(color=Progression_Status), shape=15, size=0, alpha=0) +
		guides(color = guide_legend(override.aes = list(size = 5, alpha=1))) +
        geom_boxplot(outlier.alpha = 0) + 
		guides(fill = FALSE) +
		geom_point(position=position_jitterdodge(jitter.width=0.2)) + 
		geom_signif(
			y_position = c(2, 2, 2, 2, 2, 2), 
			xmin = c(0.8, 1.8, 2.8, 4.8, 5.8, 6.8), 
			xmax = c(1.2, 2.2, 3.2, 5.2, 6.2, 7.2),
			annotation = c("*", "*", "*", "*", "*", "*"), 
			tip_length = 0) +
		ylim(-2, 2.35) +
        scale_fill_manual(values=progress_pal) +
        scale_color_manual(values=progress_pal) +
        labs(x="NLRC5 Target Genes", y = "Residual Exprsesion") +
        theme_classic() +
        theme(axis.text.x=element_text(angle = -45, hjust = 0, size = 10),
            legend.position="bottom", 
			legend.title = element_blank(),
			panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank())
dev.off()

#
# Figure 4A
# ------------------------------------------------------------------------------
index_bx <- readRDS("../output/PCGA_miRNA_index_BX.rds")
mrna <- readRDS("../output/PCGA_mrna_residual_BX.rds")
mirna <- readRDS("../output/PCGA_miRNA_Residual_BX.rds")

mrna.table <- data.frame(hgnc_symbol=c("CD3g", "CD19", "CD68", 
					"KRT5", "FOXJ1", "MUC5AC", "SCGB1A1", "MUC5B", "CEACAM5"),
EnsemblID=c("ENSG00000160654", "ENSG00000177455", "ENSG00000129226",
		"ENSG00000186081", "ENSG00000129654", "ENSG00000215182", "ENSG00000149021",
		"ENSG00000117983", "ENSG00000105388"))
row.names(mrna.table) <- mrna.table$EnsemblID

mir.list <- c("hsa-miR-34b-5p", "hsa-miR-449c-5p", "hsa-miR-150-5p", "hsa-miR-149-5p")
mrna.mat <- mrna[as.character(mrna.table$EnsemblID), ]
row.names(mrna.mat) <- mrna.table$hgnc_symbol
mir.mat <- mirna[mir.list, ]

out.mat <- merge(t(mir.mat), t(mrna.mat), by="row.names")
row.names(out.mat) <- out.mat[,1]
out.mat <- out.mat[,-1]

n <- nrow(index_bx)
c <- cor(out.mat)
c <- c[1:4, 5:13]
p <- cor2pvalue(c, n)
fdr.mat <- p.adjust(p, "fdr")
dim(fdr.mat) <- c(nrow(p), ncol(p))
row.names(c) <- gsub("hsa-", "", row.names(c))
dimnames(fdr.mat) <- dimnames(c)

pdf("../output/Figure_4/Figure_4A.pdf")
  corrplot(c, is.corr=FALSE,
          tl.col = "black", tl.srt = 45,
          p.mat = fdr.mat, sig.level = c(.001, .01, .05), pch.cex = 1,
          insig = "label_sig", pch.col = "white",
          col = rev(brewer.pal(n = 8, name = "RdBu")))
dev.off()
