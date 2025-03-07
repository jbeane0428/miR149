# Setup working environment
library(Biobase)
library(SummarizedExperiment)
library(RColorBrewer)
library(wesanderson)
library(GGally)
library(tidyverse)
library(data.table)
library(gridExtra)
library(grid)
library(fgsea)
library(ggridges)

# ==========================================================================
# Figure 4B: expression enrichment
# ==========================================================================
compartment.pal <- wes_palette("Darjeeling1")[c(1,4)]
names(compartment.pal) <- c("Epithelial", "Immune")

# Read and process FANTOM5 miRNA data
fantom.index <- as.data.frame(fread("../input/human.srna.samples.tsv"))
row.names(fantom.index) <- fantom.index$name
row.names(fantom.index) <- gsub("\\+", "\\.", row.names(fantom.index))
fantom.index <- fantom.index[, -1]
fantom.count <- read.table("../input/human.srna.cpm_include_novel.txt", header=T, sep='\t', row.names=1)

# Get Cell Ontology formatted
onto <- read.table("../input/human.srna.cellontology.tsv", header=T, sep='\t')
names(onto) <- c("Ontology", "ID")
onto$Ontology <- tolower(onto$Ontology)
onto.gs <- split(as.character(onto$ID), as.factor(onto$Ontology))
onto.gs <- lapply(onto.gs, function(x) unlist(str_split(x, ",")))

# Get miR-149-5p expression
out.mat <- merge(fantom.index, t(fantom.count["hsa-miR-149-5p", ]), by="row.names")
names(out.mat)[4] <- "miR-149-5p"
out.mat$"miR-149-5p" <- log2(out.mat$"miR-149-5p" + 0.1)

# Get cell type compartment
# Using "epithelial cell" and "leukocyte"
out.mat$Immune <- ifelse(out.mat$Row.names %in% onto.gs[["leukocyte"]], "Immune", "Others")
out.mat$Epithelial <- ifelse(out.mat$Row.names %in% onto.gs[["epithelial cell"]], "Epithelial", "Others")
out.mat$Compartment <- out.mat$Epithelial
out.mat[out.mat$Immune == "Immune", ]$Compartment <- "Immune"
table(out.mat$Compartment)

#
# Test of enrichment by compartments
# ------------------------------------------------------------------------------
gs <- list(
		out.mat[out.mat$Epithelial != "Others", ]$Row.names,
		out.mat[out.mat$Immune != "Others", ]$Row.names)
names(gs) <- c("Epithelial", "Immune")
cell.list <- setNames(out.mat[,4], out.mat$Row.names)

set.seed(1234)
fgseaRes <- fgsea(pathways = gs, 
                  stats = cell.list,
                  minSize=10,
                  maxSize=4000,
				  nperm=10000)

#
# Plot
# --------------------------------------------------------------------------
out.mat <- out.mat[order(-out.mat[,"miR-149-5p"]), ]
out.mat$Position <- seq(1:nrow(out.mat))

plot.dat <- out.mat[,c("Position", "miR-149-5p", "Epithelial", "Immune")]
plot.dat <- gather(plot.dat, Name, Value, Epithelial:Immune)
plot.dat$Value <- ifelse(plot.dat$Value == "Others", NA, plot.dat$Value)
plot.dat$Name <- factor(plot.dat$Name, levels=c("Epithelial", "Immune"))
plot.dat$Value <- factor(plot.dat$Value, levels=c("Epithelial", "Immune"))

pdf("../output/Figure_4/Figure_4B.pdf", width=5, height=5)
	p1 <- ggplot(plot.dat, aes(x=Position, y="miR-149-5p")) + 
			geom_segment(aes(x=plot.dat$Position, xend=plot.dat$Position, y=0, yend=plot.dat$"miR-149-5p"), size=1.2, color="grey50") +
			theme_bw() +
			xlab("") +
			ylab("miR-149-5p\nNormalized Expression") + 
			coord_fixed(ratio = 8) +
			theme(legend.position="none", 
				legend.title = element_blank(),
				axis.text.x = element_blank(), 
				axis.ticks.x = element_blank(),
				axis.title.y = element_text(size=8),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				plot.margin = unit(c(1,1,-0.5,1), "line"))

	p2 <- ggplot(plot.dat, aes(x=Position, y=Name, fill=Value)) +
			geom_tile() +
			scale_fill_manual(values = compartment.pal, na.translate = FALSE) +
			coord_fixed(ratio = 20) +
			xlab("Position") + ylab("Cell Types") + 
			guides(fill=guide_legend(title="")) +
			scale_y_discrete(limits=rev) +
			theme_ridges(grid = FALSE) + 
			theme_bw() +
			theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				plot.margin = unit(c(-0.5,1,1,1), "line"), 
				axis.title.y = element_text(size=8),
				axis.title.x = element_text(size=8),
				axis.text.y=element_text(size = 6),
				legend.position='bottom')

	grid.newpage()
	grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
dev.off()


# ==========================================================================
# Figure 4C: correlation by cell compartments
# ==========================================================================
out.mat$ID <- substr(out.mat$Row.names, 18, nchar(out.mat$Row.names))
out.mat$ID <- gsub("\\.", "+", out.mat$ID)

fantom.mrna <- read.table("../input/human.cage.tpm.txt", header=T, sep='\t', check.names=F, row.names=1)

nlrc5 <- fantom.mrna["p1@NLRC5",]
nlrc5 <- nlrc5[,names(nlrc5) %in% out.mat$ID]
mrna.count <- data.frame(NLRC5 = as.numeric(nlrc5), ID=names(nlrc5))

out.mat2 <- merge(out.mat, mrna.count, by="ID", all.x=T)
out.mat2 <- out.mat2[!is.na(out.mat2$NLRC5), ]
out.mat2$NLRC5 <- log2(out.mat2$NLRC5 + 0.1)

epi.dat <- out.mat2[out.mat2$Epithelial == "Epithelial", ]
cor.test(epi.dat$"miR-149-5p", epi.dat$NLRC5)
# t = -7.1289, df = 85, p-value = 3.095e-10
#        cor
# -0.6116997

imm.dat <- out.mat2[out.mat2$Immune == "Immune", ]
cor.test(imm.dat$"miR-149-5p", imm.dat$NLRC5)
# t = -1.68, df = 34, p-value = 0.1021
#        cor
# -0.2768614

pdf("../output/Figure_4/Figure_4C.pdf", width=4, height=4)
	use.label <- c(paste0("italic(R) ==", "-0.61"),
                    paste0("italic(p)-value <", "0.001"))
    ggplot(epi.dat, aes(x=`miR-149-5p`, y=NLRC5)) +
        geom_point(size=2, color=compartment.pal["Epithelial"]) +
        geom_smooth(method=lm, 
					size = 0.5, 
					linetype = "dashed",
					color=compartment.pal["Epithelial"], 
					fill=compartment.pal["Epithelial"],
					alpha = 0.25) + 
		annotate("text", x = max(epi.dat$`miR-149-5p`), y = c(max(epi.dat$NLRC5)*1.2, max(epi.dat$NLRC5)*1.08),
				label = use.label, parse = TRUE, hjust=1) +
		xlab("miR-149-5p\nNormalized Expression") + ylab("NLRC5\nNormalized Expression") +
        theme_bw() + 
		theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank())
dev.off()

