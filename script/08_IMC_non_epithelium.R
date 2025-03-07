library(tidyverse)
library(data.table)
library(RColorBrewer)
library(cowplot)
library(ggrepel)
library(ggpubr) 
library(ggsignif)
library(rstatix)
library(ComplexHeatmap)
library (lme4)
library (lmerTest)
library(emmeans) 
library(SpatialExperiment)
library(imcRtools)
library(Seurat)
library(dittoSeq)
library(patchwork)
library(sccomp)


# define cell types order and colors
StromaCell_order <- c("B cells", "CD4 T cells", "CD8 T cells", "CD8 T/CD4 T/M1 macrophage cells",
                      "M1 macrophages", "M2 macrophages", "M1/M2 macrophages",  "NK cells/Neutrophil",  
                      "SMA+ stroma cells", "Vim+ stroma cells", "SMA+Vim+ stroma cells",  "Undefined stroma")

StromaCell_color <- setNames(c( "#007ab9",      "#00b03a",    "#71e9b1",  "#0092A7",    
                                "#d93fbc",   "#eb8adf",   "#7D3C98", "#D7BDE2",  
                                "#f9bf86",  "#dcc5a6",  "#c05d27",   "gray" ),
                             c("B cells", "CD4 T cells", "CD8 T cells", "CD8 T/CD4 T/M1 macrophage cells",
                               "M1 macrophages", "M2 macrophages", "M1/M2 macrophages",  "NK cells/Neutrophil",  
                               "SMA+ stroma cells", "Vim+ stroma cells", "SMA+Vim+ stroma cells",  "Undefined stroma"))


# load spe for cells in non-epithelium
spe_stroma <- readRDS( "../input/spe_stroma.rds")

metadata_stroma <- as.data.frame(colData(spe_stroma))
metadata_stroma$CellType_Stroma <- factor(metadata_stroma$CellType_Stroma, levels = StromaCell_order)
metadata_stroma$Progression_Status = factor(metadata_stroma$Progression_Status, levels= c("Regressive", "Progressive/ \n persistent"))
colData(spe_stroma) <- DataFrame(metadata_stroma)
metadata(spe_stroma)$color_vectors$CellTypes_stroma <- StromaCell_color
# saveRDS(spe_stroma, "../input/spe_stroma.rds")

# ==========================================================================
# Figure 5F: tSNE for non-epithelial cells
# ==========================================================================
axis_plot <- ggplot(data.frame(x= 100, y = 100), aes(x = .data[["x"]], y = .data[["y"]])) +
  geom_point() +
  xlim(c(0, 10)) + ylim(c(0, 10)) +
  theme_classic() +
  ylab("tSNE2") + xlab("tSNE1") +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line( 
          arrow = arrow(angle = 15, length = unit(.5, "cm"), type = "closed"))  )

figure_layout <- c(
  area(t = 1, l = 1, b = 11, r = 11),
  area(t = 10, l = 1, b = 11, r = 2))

p <- dittoDimPlot(spe_stroma, var = "CellType_Stroma",
                  reduction.use = "TSNE_log10Zstroma", size = 0.2, do.label = FALSE,
                  show.axes.numbers = FALSE, ylab = "tSNE2", xlab = "tSNE1",
                  main = "", 
                  legend.title = "Cell Types in Stroma", legend.size = 8)+
  scale_color_manual(values = metadata(spe_stroma)$color_vectors$CellTypes_stroma, name = "Cell Types in Non-epithelium") + 
  theme_void() + axis_plot + plot_layout(design = figure_layout)

pdf("../output/Figure_5/Figure_5F.pdf", height = 5, width = 7)
print(p)
dev.off()


# ==========================================================================
# S.Figure 8A: TNSE for non-epithelium
# ==========================================================================
# visualize Sample ID
p1 <- dittoDimPlot(spe_stroma, var = "SampleID", 
                   reduction.use = "TSNE_log10Zstroma", size = 0.02, do.label = FALSE,
                   show.axes.numbers = FALSE, ylab = "tSNE2", xlab = "tSNE1",
                   main = "",
                   legend.title = "Sample ID", legend.size = 6) + # legend.title doesn't change the legend title
  scale_color_manual(values = metadata(spe_stroma)$color_vectors$SampleID, name = "Sample ID") + 
  theme_void() + axis_plot + plot_layout(design = figure_layout)

# dir.create("../output/SFigure_8")
pdf("../output/SFigure_8/SFigure_8A1.pdf", height = 4, width = 4.5)
print(p1)
dev.off()

p2 <- dittoDimPlot(spe_stroma, var = "Progression_Status", 
                   reduction.use = "TSNE_log10Zstroma", size = 0.02,
                   main = "",
                   legend.title = "Progression Status", legend.size = 6) + 
  scale_color_manual(values = metadata(spe_stroma)$color_vectors$Progression_Status, name = "Progression Status") +
  theme_void() + axis_plot + plot_layout(design = figure_layout)

pdf("../output/SFigure_8/SFigure_8A2.pdf", height = 4, width = 5)
print(p2)
dev.off()

p3 <- dittoDimPlot(spe_stroma, var = "Genomic_Smoking_Status", 
                   reduction.use = "TSNE_log10Zstroma", size = 0.02,
                   main = "",
                   legend.title = "Smoking Status", legend.size = 6) + 
  scale_color_manual(values = metadata(spe_stroma)$color_vectors$Smoking_Status, name = "Smoking Status") +
  theme_void() + axis_plot + plot_layout(design = figure_layout)

pdf("../output/SFigure_8/SFigure_8A3.pdf", height = 4, width = 4.5)
print(p3)
dev.off()


# ==========================================================================
# Figure 5G: Heatmap of average expression of cell types in non-epithelium
# ==========================================================================
# compute the log-Z normalized within non-epithelium count matrix
count_matrix <- as.data.frame(t(assay(spe_stroma, "log10")))
metadata_cell <- as.data.frame(colData(spe_stroma))
count_matrix$ObjectID = rownames(count_matrix)
count_matrix$image_id = metadata_cell$image_id
image_ids = unique(count_matrix$image_id)

countMatrix_Log_Z = data.frame()
subset_data = c("image_id", "ObjectID")

for (id in image_ids){
  countMatrix_temp <- count_matrix %>% filter(image_id == id)
  celldata_subset <- countMatrix_temp[,subset_data]
  countMatrix_temp = countMatrix_temp[,!names(countMatrix_temp) %in% subset_data]
  # Call z-score normalize
  countMatrix_temp <- scale(countMatrix_temp,center=T,scale=T)
  countMatrix_temp = cbind(countMatrix_temp, celldata_subset)
  countMatrix_Log_Z<- rbind(countMatrix_Log_Z, countMatrix_temp)
}

all.equal(countMatrix_Log_Z$ObjectID, rownames(metadata_cell))
countMatrix_Log_Z1 = countMatrix_Log_Z[,!names(countMatrix_Log_Z) %in% subset_data]
assay(spe_stroma, "log10_Z_stroma") <- t(countMatrix_Log_Z1)

# saveRDS(spe_stroma, "../input/spe_stroma.rds")

# get the count matrix for heatmap
expMatDF_stroma <- as.data.frame(t(assay(spe_stroma, "log10_Z_stroma")))
expMatDF_stroma_logZstromalNormalized <- as.data.frame(t(assay(spe_stroma, "log10_Z_stromalNormalized")))

# Set gene list to be displayed on heatmap
Canonical_stroma_markers_for_heatmap = 
  c("CD20", "CD3", "CD4","CD8a",  "GZMB", "CD66b", "CD68", "CD163",   "HLA_DR", "PDL1",
    "aSMA", "Vimentin") 

NLRC5_target_markers_for_heatmap = c("NLRC5", "B2M",  "PSMB9", "TAP1") 

# Add cluster information to data frame
expMatDF_stroma$Cluster = spe_stroma$CellType_Stroma
expMatDF_stroma_logZstromalNormalized$Cluster = spe_stroma$CellType_Stroma

# Take average expression of feature per cluster
# canonical markers
expMatDF_stroma_SummaryMean <- expMatDF_stroma %>%
  group_by(Cluster) %>%
  summarise(across(
    .cols = all_of(Canonical_stroma_markers_for_heatmap),
    .fns = list(
      mean = ~ mean(.x, na.rm = TRUE)))) %>% as.data.frame()
rownames(expMatDF_stroma_SummaryMean) <- expMatDF_stroma_SummaryMean$Cluster
expMatDF_stroma_SummaryMean$Cluster <- NULL


# NLRC5 markers
expMatDFSummaryMean_stroma_logZstromalNormalized <- expMatDF_stroma_logZstromalNormalized %>%
  group_by(Cluster) %>%
  summarise(across(
    .cols = all_of(NLRC5_target_markers_for_heatmap),
    .fns = list(
      mean = ~ mean(.x, na.rm = TRUE)))) %>% as.data.frame()
rownames(expMatDFSummaryMean_stroma_logZstromalNormalized) <- expMatDFSummaryMean_stroma_logZstromalNormalized$Cluster
expMatDFSummaryMean_stroma_logZstromalNormalized$Cluster <- NULL

# re-order
expMatDF_stroma_SummaryMean = expMatDF_stroma_SummaryMean[StromaCell_order,]
expMatDFSummaryMean_stroma_logZstromalNormalized = expMatDFSummaryMean_stroma_logZstromalNormalized[StromaCell_order,]

colnames(expMatDF_stroma_SummaryMean) <- gsub("_mean", "", colnames(expMatDF_stroma_SummaryMean))
colnames(expMatDFSummaryMean_stroma_logZstromalNormalized) <- gsub("_mean", "", colnames(expMatDFSummaryMean_stroma_logZstromalNormalized))

# cell type proportion
propTableForHeatmap_stroma <- as.matrix(table(expMatDF_stroma$Cluster)) %>% prop.table(2) *100
propTableForHeatmap_stroma <- propTableForHeatmap_stroma[StromaCell_order,]

# heatmap annotation
rowAnno_right_stroma = ComplexHeatmap::rowAnnotation(
  `% Cell Types`  = ComplexHeatmap::anno_barplot(
     propTableForHeatmap_stroma, which = "row", axis = TRUE, bar_width = 0.8, width = unit(2.5, "cm"),
     gp = grid::gpar(fill = StromaCell_color), 
     border = TRUE,
     show_annotation_name = FALSE))

rowAnno_left_stroma = ComplexHeatmap::rowAnnotation(
  `Cell Types` = StromaCell_order, 
  col = list(`Cell Types` = StromaCell_color),
  simple_anno_size = unit(0.3, "cm"),
  show_legend = FALSE,
  show_annotation_name = TRUE,
  annotation_name_rot = 90
)

# Create color scale for heatmap
col1 = circlize::colorRamp2(c(-3, 0, 3),
                            c("#1E90FF", "#FFFFFF", "#CD2626"))

col2 = circlize::colorRamp2(c(-1, 0, 2),
                            c("#02f900", "#FFFFFF", "#ffa921"))

# Create heatmap
ht1_stroma <- ComplexHeatmap::Heatmap(matrix = expMatDF_stroma_SummaryMean, 
                                      heatmap_legend_param = list(
                                        col_fun = col1,
                                        title = "Average Expression\n(Log-Z normalized within non-epithelium)", 
                                        direction = "horizontal"),
                                      col = col1,
                                      cluster_rows = FALSE,
                                      cluster_columns = FALSE, column_names_rot = 90,
                                      left_annotation = rowAnno_left_stroma, 
                                      row_names_gp = grid::gpar(fontsize = 10), 
                                      column_names_gp = grid::gpar(fontsize = 10), 
                                      show_column_names = TRUE,
                                      show_row_names = FALSE,
                                      column_title = "Canonical phenotypic markers"
)

ht2_stroma <- ComplexHeatmap::Heatmap(matrix = as.matrix(expMatDFSummaryMean_stroma_logZstromalNormalized),
                                      heatmap_legend_param = list(
                                        col_fun = col2,
                                        title = "Average Expression\n(Log-Z stromal normalized)",
                                        direction = "horizontal"),
                                      col= col2,
                                      cluster_rows = FALSE,
                                      cluster_columns = FALSE, column_names_rot = 90,
                                      right_annotation = rowAnno_right_stroma,
                                      row_names_gp = grid::gpar(fontsize = 10), 
                                      column_names_gp = grid::gpar(fontsize = 10), 
                                      show_column_names = TRUE,
                                      show_row_names = FALSE,
                                      column_title = "NLRC5 targets")

p<- ComplexHeatmap::draw(ht1_stroma + ht2_stroma , heatmap_legend_side = "bottom") 


pdf("../output/Figure_5/Figure_5G.pdf", height = 4, width = 6)
print(p)
dev.off()


# ==========================================================================
# S.Figure 8B: stack plot of sample ID across cell types in non-epithelium
# ==========================================================================
cells_df <- colData(spe_stroma) %>% as.data.frame() %>%
  mutate(Progression_Status = factor(Progression_Status, levels= c("Regressive", "Progressive/ \n persistent")))

# stackplot : y = sample ID percentage, x = cell types
cells_table_df <- table(cells_df$CellType_Stroma, cells_df$SampleID) %>% as.data.frame()

p <- ggplot(cells_table_df, aes(y = Freq)) +
  labs(x = "Cell Types in Stroma", y = "Proportion (%)", title = "") + 
  geom_bar( aes(x = Var1, fill = Var2), 
            position = "fill", stat = "identity") + # 
  scale_fill_manual(values = metadata(spe_epi)$color_vectors$SampleID, name = "Sample ID") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
  theme_bw() + theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, angle = -45, hjust = 0, vjust = 1),
    axis.title = element_text(size = 14, color = "black"),
    legend.key.height  =  unit(0.8, "lines"),
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 12),
    panel.border = element_blank(),
    panel.spacing.x = unit(1, "lines")) 


pdf("../output/SFigure_8/SFigure_8B.pdf", height = 5, width = 6)
print(p)
dev.off()


# ==========================================================================
# Figure 5H: Differential NLRC5 expression between cells in epithleium and non-epithleium
# ==========================================================================
spe_all <- readRDS( "../input/spe_all.rds")

count_matrix_all <- as.data.frame(t(assay(spe_all, "log10_Z_stromalNormalized")))
metadata_all <- as.data.frame(colData(spe_all))
all.equal(rownames(metadata_all), rownames(count_matrix_all))

metadata_all$NLRC5_logZ_stromalNormalized = count_matrix_all$NLRC5
metadata_all$TAP1_logZ_stromalNormalized = count_matrix_all$TAP1
metadata_all$HLA_ABC_logZ_stromalNormalized = count_matrix_all$HLA_ABC
metadata_all$PSMB9_logZ_stromalNormalized = count_matrix_all$PSMB9
metadata_all$B2M_logZ_stromalNormalized = count_matrix_all$B2M
colData(spe_all) <- DataFrame(metadata_all)

p= ggplot(metadata_all,aes( x = isEpithelium_byCluster, y = NLRC5_logZ_stromalNormalized, fill= isEpithelium_byCluster)) + 
  labs(x = NULL, y = "Log-Z stromal normalized NLRC5",
       title = "") + 
  theme_classic() +  
  geom_violin(width=1, show.legend = FALSE) +
  geom_boxplot(width=0.3, color="darkgrey", alpha=1, show.legend = FALSE, outlier.shape = NA) + 
  scale_fill_manual(values = c("red", "orange")) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, size=12, color = "black"),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  stat_compare_means(comparisons = list( c("Epithelium", "Non_epithelium") ), label = "p.signif", method = "wilcox.test") # Add pairwise comparisons p-value

print(p)

pdf("../output/Figure_5/Figure_5H.pdf", height = 4, width = 3)
print(p)
dev.off()


# ==========================================================================
# S.Figure 8C: stack plot of cell types between progression status
# ==========================================================================
cells_df <- colData(spe_stroma) %>% as.data.frame() %>%
  mutate(Progression_Status = factor(Progression_Status, levels= c("Regressive", "Progressive/ \n persistent")))

cells_table_df <- table(cells_df$SampleID, cells_df$CellType_Stroma) %>% as.data.frame() 
colnames(cells_table_df) <- c("Sample", "CellType",  "Freq")

cells_table_df$Progression_Status <- cells_df$Progression_Status[match(cells_table_df$Sample, cells_df$SampleID)]

p <- ggplot(cells_table_df, aes(y = Freq)) +
  labs(x = "Sample ID", y = "Proportion (%)", title = "") + 
  facet_wrap("Progression_Status", scales = "free_x") +
  geom_bar( aes(x = Sample, fill = CellType), 
            position = "fill", stat = "identity") + # 
  scale_fill_manual(values = StromaCell_color, name = "Cell types in non-epithelium") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, face = "bold"), # facet title
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 14, color = "black"),
    legend.key.height  =  unit(0.8, "lines"),
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 12),
    panel.border = element_blank(),
    panel.spacing.x = unit(1, "lines")) 


pdf("../output/SFigure_8/SFigure_8C.pdf", height = 5, width = 8.5)
print(p)
dev.off()


# ==========================================================================
# Figure 6B: Differential cell type composition between progression status
# ==========================================================================
# run sccomp
res_stroma_fComp_prog =
  spe_stroma |>
  sccomp_estimate(
    formula_composition = ~ Progression_Status,
    .sample = SampleID,
    .cell_group = CellType_Stroma, 
    mcmc_seed = 12345,
    variational_inference = FALSE,
    cores = 1) |>
  sccomp_remove_outliers(variational_inference = FALSE, mcmc_seed = 12345)  |>
  sccomp_test(test_composition_above_logit_fold_change = 0.1)


# saveRDS(res_stroma_fComp_prog,  "../output/sccomp_CellType_Stroma_fComposition_Progression.rds")
# res_stroma_fComp_prog <- readRDS(  "../output/sccomp_CellType_Stroma_fComposition_Progression.rds")

res_stroma_fComp_prog1 <- res_stroma_fComp_prog %>% filter(parameter != "(Intercept)") %>% 
  select(CellType_Stroma, parameter, c_effect, c_pH0, c_FDR) %>%
  mutate(c_FDR_crop = ifelse(c_FDR< 0.0001, 0.0001, c_FDR)) 

vplot <- ggplot(res_stroma_fComp_prog1) +
  aes(y=-log10(c_FDR_crop), x=c_effect, label =  CellType_Stroma) +
  scale_color_manual(values = StromaCell_color)+
  geom_point(size=5, aes(color = CellType_Stroma), show.legend = FALSE) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", color="red", linewidth=0.5) +
  annotate("text", x=0.5, y=-log10(0.05)+0.1,
           label=paste("FDR<0.05"), size=5, fontface="bold", color = "red") +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", linewidth=0.5) +
  geom_text_repel(size = 4) +
  labs(y="-Log10(FDR)", x = "Estimate, Progressive - Regressive") + 
  theme_bw()+
  theme(axis.text.x = element_text( size=16, color = "black"),
        axis.text.y = element_text(size=16, color = "black"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18)
  ) 


pdf("../output/Figure_6/Figure_6B.pdf.pdf", height = 5, width = 5.7)
print(vplot)
dev.off()




