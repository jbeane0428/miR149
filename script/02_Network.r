# Setup working environment
library(Biobase)
library(SummarizedExperiment)
library(biomaRt)
library(genefilter)
library(ggplot2)
library(gplots)
library(limma)
library(RColorBrewer)
library(sva)
library(GGally)
library(network)
library(sna)
library(GSVA)
library(tidyverse)
library(ggridges)
library(reshape2)
library(corrplot)
library(igraph)
library(ggridges)
library(qgraph)

# ==============================================================================
# Necessary functions
# ==========================================================================
cor2pvalue = function(r, n) {
	t <- (r*sqrt(n-2))/sqrt(1-r^2)
	p <- 2*(1 - pt(abs(t),(n-2)))
	return(p)
}

Run.module.f.test <- function(target.matrix, binary.cor.matrix, module.gene, connect.mirna) {
    target.module <- target.matrix[row.names(target.matrix) %in% module.gene,]
    target.others <- target.matrix[!row.names(target.matrix) %in% module.gene,]

    pval.cov.bin.module <- binary.cor.matrix[row.names(binary.cor.matrix) %in% module.gene,]
    pval.cov.bin.others <- binary.cor.matrix[!row.names(binary.cor.matrix) %in% module.gene,]

	print(all.equal(dimnames(pval.cov.bin.module), dimnames(target.module)))
	print(all.equal(dimnames(pval.cov.bin.others), dimnames(target.others)))
    target.cor.module <- pval.cov.bin.module + target.module
    target.cor.others <- pval.cov.bin.others + target.others

	fisher.table <- matrix(,nrow=0, ncol=3)
    for (mir in connect.mirna) {
		out <- c()

        this.module <- target.cor.module[,mir]
        this.others <- target.cor.others[,mir]

        this.module.target <- target.module[,mir]
        this.others.target <- target.others[,mir]
    
        a <- sum(this.module==2)
        b <- sum(this.others==2)
        c <- sum(this.module.target==1) - a
        d <- sum(this.others.target==1) - b

        contingency.table <- matrix(c(a,c,b,d), nrow = 2)
        res <- fisher.test(contingency.table)

        p <- res$p.value
        or <- res$estimate
		
        out <- c(mir, or, p)

        fisher.table <- rbind(fisher.table, out)
    }

    fisher.table <- as.data.frame(fisher.table)
    fisher.table[,2] <- as.numeric(as.character(fisher.table[,2]))
    fisher.table[,3] <- as.numeric(as.character(fisher.table[,3]))
    fisher.table$FDR <- p.adjust(fisher.table[,3], "fdr")
    row.names(fisher.table) <- fisher.table[,1]
    fisher.table <- fisher.table[,-1]
	
	return(fisher.table)
}

Run.module.t.test <- function(target.matrix, correlation.matrix, module.gene, connect.mirna) {
    target.module <- target.others <- target.matrix
    target.module[!row.names(target.module) %in% module.gene,] <- 0
    target.others[row.names(target.others) %in% module.gene,] <- 0

    # T test comparing mean correlation difference
	print(all.equal(dimnames(correlation.matrix), dimnames(target.module)))
	print(all.equal(dimnames(correlation.matrix), dimnames(target.others)))
    target.cor.mat <- t(correlation.matrix * target.module)
    nontarget.cor.mat <- t(correlation.matrix * target.others)
    target.cor.mat[target.cor.mat == 0] <- NA
    nontarget.cor.mat[nontarget.cor.mat == 0] <- NA

	# Fisher Z transformation
    target.cor.mat <- log((1+target.cor.mat)/(1-target.cor.mat))/2
    nontarget.cor.mat <- log((1+nontarget.cor.mat)/(1-nontarget.cor.mat))/2

	select.mir <- c()
    t.table <- matrix(,nrow=0, ncol=3)
    for (mir in connect.mirna) {
        temp.t <- target.cor.mat[mir, ]
        d.t <- data.frame(Z = temp.t, Group = "Target")
        d.t <- d.t[!is.na(d.t$Z), ]
        temp.n <- nontarget.cor.mat[mir, ]
        d.n <- data.frame(Z = temp.n, Group = "Non Target")
        d.n <- d.n[!is.na(d.n$Z), ]
        
		out <- c()
		# Making the t-stat value negative for left shifting of in-module targets PCC
       if (nrow(d.t) > 1 & nrow(d.n) > 1 ) {
            out.table <- rbind(d.t, d.n)
			out.table$Group <- factor(out.table$Group, levels=c("Target", "Non Target"))
            res <- t.test(out.table$Z ~ out.table$Group)
            out <- c(mir, res$statistic, res$p.value)
            t.table <- rbind(t.table, out)
        } else if (nrow(d.t) == 1 & nrow(d.n) > 1) {
            res <- t.test(d.n$Z, mu=d.t$Z)
            out <- c(mir, -res$statistic, res$p.value)
            t.table <- rbind(t.table, out)
        } else if (nrow(d.t) > 1 & nrow(d.n) == 1) {
            res <- t.test(d.t$Z, mu=d.n$Z)
            out <- c(mir, res$statistic, res$p.value)
            t.table <- rbind(t.table, out)
        } else { # No miRNA had 0 targets in other module but >=1 targets in module of interests
            out <- c(mir, NA, NA)
            t.table <- rbind(t.table, out)
        }
    }
    t.table <- as.data.frame(t.table)
    t.table[,2] <- as.numeric(as.character(t.table[,2]))
    t.table[,3] <- as.numeric(as.character(t.table[,3]))
    t.table$FDR <- p.adjust(t.table[,3], "fdr")
    row.names(t.table) <- t.table[,1]
    t.table <- t.table[,-1]
	
	return(t.table)
}

# ==============================================================================
# Calculate correlation matrix
# ==============================================================================
mrna <- readRDS("../output/PCGA_mrna_residual_BX.rds")
mirna <- readRDS("../output/PCGA_miRNA_Residual_BX.rds")

mrna <- mrna[, colnames(mirna)]

mat <- rbind(mrna, mirna)

n <- ncol(mat)
cov.mat <- cor(t(mat))
r <- cov.mat[1:16653, 16654:ncol(cov.mat)]
p <- cor2pvalue(r, n)

cormat.list <- list(r, p)
saveRDS(cormat.list, "../output/BX_Residual_correlation_matrix.rds")

# ==============================================================================
# Construct miRNA-gene network and filter
# ==============================================================================
gene_module <- readRDS("../input/PCGA_Gene_Module.rds")
g_list <- row.names(gene_module)

# Preprocess target and correlation matrix to include only module genes
target.combine <- readRDS("../input/Combine_Target_Matrix.rds")
target.combine <- target.combine[row.names(target.combine) %in% g_list, ]

cov.mat <- cormat.list[[1]]
pval.mat <- cormat.list[[2]]
cov.mat <- cov.mat[row.names(cov.mat) %in% g_list, ]
pval.mat <- pval.mat[row.names(pval.mat) %in% g_list, ]

# Filter by correlation strength and significance
# Generate matrix with significant negative correlation
cov.mat.bin <- ifelse(cov.mat < 0, 1, 0)
fdr.mat <- pval.mat
fdr.mat[] <- p.adjust(pval.mat, "fdr")
pval.mat.bin <- ifelse(fdr.mat < 0.05, 1, 0)
pval.cov.mat.bin <- cov.mat.bin + pval.mat.bin
pval.cov.mat.bin <- ifelse(pval.cov.mat.bin == 2, 1, 0)

#
# Generate network significant negative correlation between miRNAs and predicted targets
# ------------------------------------------------------------------------------
all.equal(dimnames(target.combine), dimnames(pval.cov.mat.bin))
cov.sum.mat <- t(target.combine) + t(pval.cov.mat.bin)
cov.sum.mat <- ifelse(cov.sum.mat == 2, 1, 0)
cov.sum.mat <- cov.sum.mat[rowSums(cov.sum.mat == 1) != 0, colSums(cov.sum.mat == 1) != 0]

#
# Generate full network between miRNA and modules
# If any significant negative target-pair existed, the connection is retained
# Between the miRNA and the corresponded module
# --------------------------------------------------------------------------
net.mat <- cov.sum.mat
net.mat.module <- matrix(nrow=nrow(net.mat), ncol=9)
colnames(net.mat.module) <- names(table(gene_module$Gene_Module_Number))
row.names(net.mat.module) <- row.names(net.mat)

for (i in names(table(gene_module$Gene_Module_Number))) {
    net.mat.module[,i] <- ifelse(rowSums(net.mat[,colnames(net.mat) %in% gene_module[gene_module$Gene_Module_Number == i,]$EnsemblID]) > 0, 1, 0)
}

#
# Filtering the connection
# The major analysis: applying density shift and target enrichment analysis
# to retain miRNA-to-module connection that is more specific for each module
# ------------------------------------------------------------------------------
sig.cutoff <- 0.05
mir.module.gs <- t.res.list <- f.res.list <- list()

for (Module in names(table(gene_module$Gene_Module_Number))) {
	# First identified the list of miRNAs with significant negative connection to target in the module
    mir_list <- row.names(net.mat.module)[net.mat.module[, Module] == 1]

	# ==============================================================================
	# F test
	# ==============================================================================
	fisher.res <- c()
	fisher.res <- Run.module.f.test(target.matrix = target.combine,
								binary.cor.matrix = pval.cov.mat.bin,
								module.gene = row.names(gene_module[gene_module$Gene_Module_Number == Module,]),
								connect.mirna = mir_list)
	f.res.list[[Module]] <- fisher.res
	f.mir <- row.names(fisher.res[fisher.res$FDR < sig.cutoff & fisher.res[,1] > 1, ])

	# ==============================================================================
	# T
	# ==============================================================================
	t.res <- c()
	t.res <- Run.module.t.test(target.matrix = target.combine,
								correlation.matrix = cov.mat,
								module.gene = row.names(gene_module[gene_module$Gene_Module_Number == Module,]),
								connect.mirna = mir_list)
	t.res.list[[Module]] <- t.res
	t.mir <- row.names(t.res[t.res$FDR < sig.cutoff & t.res[,1] < 0,])

	select.mir <- t.mir[t.mir %in% f.mir]
	mir.module.gs[[Module]] <- select.mir
}

saveRDS(mir.module.gs, "../output/miRNA_by_PCGA_Gene_Module.rds")

print(t.res.list[["Module9"]][mir.module.gs[["Module9"]], ])
print(f.res.list[["Module9"]][mir.module.gs[["Module9"]], ])

#
# Filter the miRNA-gene network to keep only miRNA assigned to each module
# --------------------------------------------------------------------------
net.mat.fil <- matrix(0, nrow=nrow(net.mat), ncol=ncol(net.mat))
colnames(net.mat.fil) <- colnames(net.mat)
row.names(net.mat.fil) <- row.names(net.mat)

for (Module in names(table(gene_module$Gene_Module_Number))) {

    temp.mat <- net.mat[, colnames(net.mat) %in% as.character(gene_module[gene_module$Gene_Module_Number == Module,]$EnsemblID)]
    temp.mat <- temp.mat[rowSums(temp.mat == 1)!=0, colSums(temp.mat == 1)!=0]
   
    edge.mat <- temp.mat[mir.module.gs[[Module]], ,drop=F]
    edge.list <- reshape2::melt(edge.mat)
    edge.list <- as.matrix(edge.list[edge.list$value == 1, ,drop=F])

    net.mat.fil[edge.list[,c(1,2)]] <- 1

	print(Module)
}

saveRDS(net.mat.fil, "../output/PCGA_BX_Filtered_miRNA_Gene_Network.rds")

#
#Filter the miRNA-module network to keep only miRNA assigned to each module
# --------------------------------------------------------------------------
net.mat.fil.module <- matrix(nrow=nrow(net.mat.fil), ncol=9)
colnames(net.mat.fil.module) <- names(table(gene_module$Gene_Module_Number))
row.names(net.mat.fil.module) <- row.names(net.mat.fil)

for (Module in names(table(gene_module$Gene_Module_Number))) {
    net.mat.fil.module[,Module] <- ifelse(rowSums(net.mat.fil[,colnames(net.mat.fil) %in% gene_module[gene_module$Gene_Module_Number == Module,]$EnsemblID]) > 0, 1, 0)
}

saveRDS(net.mat.fil.module, "../output/PCGA_BX_Filtered_miRNA_GeneModule_Network.rds")

# ==============================================================================
# Plotting
# ==============================================================================

#
# Figure 2A
# ------------------------------------------------------------------------------
getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
net.pal <- getPalette(9)
net.pal <- c("grey60", net.pal)
names(net.pal) <- c("miRNA", colnames(net.mat.fil.module))

net.mat.fil.module <- readRDS("../output/PCGA_BX_Filtered_miRNA_GeneModule_Network.rds")
module.graph <- net.mat.fil.module %>%
				igraph::graph_from_incidence_matrix() %>%
				igraph::delete.vertices(degree(.)==0)

V(module.graph)$type <- ifelse(V(module.graph)$name %in% colnames(net.mat.fil.module), V(module.graph)$name, "miRNA")
V(module.graph)$color <- net.pal[V(module.graph)$type]

e <- get.edgelist(module.graph,names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(module.graph))
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

V(module.graph)$name <- ifelse(V(module.graph)$name %in% colnames(net.mat.fil.module), V(module.graph)$name, NA)

pdf("../output/Figure_2/Figure_2A.pdf", width=6, height=6)
	par(mar=c(0,0,0,0)+.1)
	plot(module.graph,
		layout=l, 
		rescale=FALSE, 
		vertex.size=ifelse(!is.na(V(module.graph)$name), 12, 3),
		vertex.label.family="Helvetica", 
		vertex.label.cex=0.4,
		vertex.frame.color="gray25",
		vertex.label.color="black", 
		edge.color="gray80")

	legend("topleft", legend=c("Gene Module", "miRNA"), pch=21, pt.cex=c(2,1), cex=0.8, pt.bg="grey60")
dev.off()

#
# Figure 2B
# ------------------------------------------------------------------------------
mirna_res <- readRDS("../output/PCGA_miRNA_Residual_BX.rds")
modules <- readRDS("../output/PCGA_mRNA_modules_score_BX.rds")
mir.module.gs <- readRDS("../output/miRNA_by_PCGA_Gene_Module.rds")

mir.module.dat <- data.frame(miRNA = unlist(mir.module.gs), 
                            Gene_Module_Number = rep(names(mir.module.gs), time=lapply(mir.module.gs, length)))

mir.module.mat <- cbind(modules, t(mirna_res))
c <- cor(mir.module.mat)
c <- as.data.frame(c[1:9, 10:ncol(mir.module.mat)])

mir.module.dat$PCC <- c[as.matrix(mir.module.dat)[, c(2,1)]]

pdf("../output/Figure_2/Figure_2B.pdf", height=4, width=4)
    ggplot(mir.module.dat, aes(x = PCC, y = Gene_Module_Number, fill = Gene_Module_Number)) +
		geom_vline(xintercept=0, linetype="dashed", color="grey75") +
		geom_density_ridges(alpha = .9, color="white") +
		scale_fill_manual(values = net.pal) +
		scale_y_discrete(limits=rev) +
		xlab("Pearson Correlation") + ylab("Gene Module") +
		theme_bw() + 
		theme(legend.position='none',
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank())
dev.off()


#
# Figure 2C
# ------------------------------------------------------------------------------
mir.list <- gsub("hsa-", "", mir.module.gs[["Module9"]])

fil.net <- readRDS("../output/PCGA_BX_Filtered_miRNA_Gene_Network.rds")
row.names(fil.net) <- gsub("hsa-", "", row.names(fil.net))
plot.net <- fil.net[mir.list, colnames(fil.net) %in% row.names(gene_module[gene_module$Gene_Module_Number == "Module9", ]), drop=F]

# Tranlate gene names
colnames(plot.net) <- gene_module[colnames(plot.net),]$hgnc_symbol

module.graph <- plot.net %>%
				igraph::graph_from_incidence_matrix() %>%
				igraph::delete.vertices(degree(.)==0)

V(module.graph)$type <- ifelse(V(module.graph)$name %in% colnames(plot.net), "Module9", "miRNA")
V(module.graph)$color <- net.pal[V(module.graph)$type]

e <- get.edgelist(module.graph,names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(module.graph))
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

pdf("../output/Figure_2/Figure_2C.pdf", width=6, height=6)
	par(mar=c(0,0,0,0)+.1)
	plot(module.graph,
		layout=l, 
		rescale=FALSE, 
		vertex.size=ifelse(V(module.graph)$name %in% colnames(plot.net), 15, 30),
		vertex.label.family="Helvetica", 
		vertex.label.cex=ifelse(V(module.graph)$name %in% colnames(plot.net), 0.6, 0.8),
		vertex.frame.color="gray25",
		vertex.label.color=ifelse(V(module.graph)$name %in% colnames(plot.net), "white", "black"), 
		edge.color="gray80")

	legend("topleft", legend=c("Module 9 Gene", "miRNA"), pch=21, pt.cex=c(1,2), cex=0.8, pt.bg=net.pal[c(10,1)])
dev.off()