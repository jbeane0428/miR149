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



#define cell types order and colors
EpiCell_Order = c(  
  "Basal cell",              
  "Basal/differentiating cell",
  "KRT5 low basal cell",
  "MUC5AC+ goblet cell",
  "MUC5B+ secretory cell",
  "MUC5B low secretory cell",
  "CEACAM5+ peri-goblet cell",                      
  "CEACAM5+SCGB1A1+ secretory cell",                                    
  "FOXJ1+ ciliated cell")

EpiCell_Color <- setNames(c("#DE3163",     "#7D0112",     "pink4",         "#ff7400",      "#ffbb5f",    "#faf39b",  # , 
                            "#d0ea56",                  "#0092A7",        "#dfb45a"),
                          c( "Basal cell", "Basal/differentiating cell", "KRT5 low basal cell", "MUC5AC+ goblet cell",  "MUC5B+ secretory cell", "MUC5B low secretory cell",
                             "CEACAM5+ peri-goblet cell",    "CEACAM5+SCGB1A1+ secretory cell",   "FOXJ1+ ciliated cell" ))


CellInEpi_Order <- c(   "Basal cell",              
                        "Basal/differentiating cell",
                        "KRT5 low basal cell",
                        "MUC5AC+ goblet cell",
                        "MUC5B+ secretory cell",
                        "MUC5B low secretory cell",
                        "CEACAM5+ peri-goblet cell",                      
                        "CEACAM5+SCGB1A1+ secretory cell",                                    
                        "FOXJ1+ ciliated cell", 
                        "CD4 T cell in epithelium",                           
                        "CD8 T cell in epithelium",                      
                        "M1/M2 macrophage in epithelium",                                  
                        "NK cell in epithelium")

CellInEpi_Color <- setNames(c("#DE3163",     "#7D0112",     "pink4",         "#ff7400",      "#ffbb5f",    "#faf39b",  # , 
                              "#d0ea56",                  "#0092A7",        "#dfb45a",     
                              "#00b03a",                 "#71e9b1",           "#7D3C98",           "#D7BDE2", "gray" ),
                            c(   "Basal cell",  "Basal/differentiating cell",  "KRT5 low basal cell", "MUC5AC+ goblet cell", "MUC5B+ secretory cell",
                                 "MUC5B low secretory cell",
                                 "CEACAM5+ peri-goblet cell",                      
                                 "CEACAM5+SCGB1A1+ secretory cell",                                    
                                 "FOXJ1+ ciliated cell", 
                                 "CD4 T cell in epithelium",                           
                                 "CD8 T cell in epithelium",                      
                                 "M1/M2 macrophage in epithelium",                                  
                                 "NK cell in epithelium", "Undefined"))

immuneCell_inEpi = c("CD4 T cell in epithelium",
                     "CD8 T cell in epithelium", 
                     "M1/M2 macrophage in epithelium", 
                     "NK cell in epithelium")

# define color vectors for histology
Histology_cell_color <- setNames(c("grey", "#F4AD31", 
                                   "#3c68b0", 
                                   "#BF0A3D", 
                                   "#3F1B03"),
                                 c("Non-annotated", "Normal epithelium", 
                                   "Hyper-metaplasia", 
                                   "Dysplasia", 
                                   "CIS"))

Histology_cell_detail_color <- setNames(c("grey", "#F4AD31", 
                                          "#0c1cb0", "#3c68b0", 
                                          "#4ab069", "#0fd412", 
                                          "#c46a84", "#BF0A3D", "red",
                                          "#3F1B03"),
                                        c("Non-annotated", "Normal epithelium", 
                                          "Goblet cell hyperplasia", "Reserve cell hyperplasia", 
                                          "Metaplasia", "Immature squamous metaplasia", 
                                          "Mild dysplasia", "Moderate dysplasia", "Severe dysplasia", 
                                          "CIS"))

# define color vector for progression status
Progression_Status_color = c("#D3DDDC", "#446455")


# ==========================================================================
# Figure 5A: TNSE for epithelium
# ==========================================================================
# load spe
spe_epi <- readRDS("../input/spe_epi.rds")

metadata <- as.data.frame(colData(spe_epi))

# cell type annotation
metadata <- metadata %>% 
  mutate(
    CellType_Epi = case_when(
      pg_clusters_TSNElogZepi_k85 == "1" ~ "CEACAM5+SCGB1A1+ secretory cell",
      pg_clusters_TSNElogZepi_k85 == "2" ~ "MUC5B+ secretory cell",
      pg_clusters_TSNElogZepi_k85 == "3" ~ "CEACAM5+SCGB1A1+ secretory cell",
      pg_clusters_TSNElogZepi_k85 == "4" ~ "CEACAM5+ peri-goblet cell",
      pg_clusters_TSNElogZepi_k85 == "5" ~ "M1/M2 macrophage in epithelium",
      pg_clusters_TSNElogZepi_k85 == "6" ~ "KRT5 low basal cell",
      pg_clusters_TSNElogZepi_k85 == "7" ~ "CD8 T cell in epithelium",
      pg_clusters_TSNElogZepi_k85 == "8" ~ "KRT5 low basal cell",
      pg_clusters_TSNElogZepi_k85 == "9" ~ "MUC5AC+ goblet cell",
      pg_clusters_TSNElogZepi_k85 == "10" ~ "CD8 T cell in epithelium",
      pg_clusters_TSNElogZepi_k85 == "11" ~ "CD4 T cell in epithelium",
      pg_clusters_TSNElogZepi_k85 == "12" ~ "Basal cell",
      pg_clusters_TSNElogZepi_k85 == "13" ~ "Basal cell",
      pg_clusters_TSNElogZepi_k85 == "14" ~ "Basal cell",
      pg_clusters_TSNElogZepi_k85 == "15" ~ "CEACAM5+ peri-goblet cell",
      pg_clusters_TSNElogZepi_k85 == "16" ~ "Basal cell",
      pg_clusters_TSNElogZepi_k85 == "17" ~ "Basal/differentiating cell",
      pg_clusters_TSNElogZepi_k85 == "18" ~ "MUC5AC+ goblet cell",
      pg_clusters_TSNElogZepi_k85 == "19" ~ "Basal cell",
      pg_clusters_TSNElogZepi_k85 == "20" ~ "MUC5B low secretory cell",
      pg_clusters_TSNElogZepi_k85 == "21" ~ "Basal cell",
      pg_clusters_TSNElogZepi_k85 == "22" ~ "MUC5AC+ goblet cell",
      pg_clusters_TSNElogZepi_k85 == "23" ~ "Basal cell",
      pg_clusters_TSNElogZepi_k85 == "24" ~ "Basal cell",
      pg_clusters_TSNElogZepi_k85 == "25" ~ "Basal cell",
      pg_clusters_TSNElogZepi_k85 == "26" ~ "Basal cell", 
      pg_clusters_TSNElogZepi_k85 == "27" ~ "CD8 T cell in epithelium",
      pg_clusters_TSNElogZepi_k85 == "28" ~ "MUC5B+ secretory cell",
      pg_clusters_TSNElogZepi_k85 == "29" ~ "Basal/differentiating cell",
      pg_clusters_TSNElogZepi_k85 == "30" ~ "Basal cell",
      pg_clusters_TSNElogZepi_k85 == "31" ~ "MUC5AC+ goblet cell",
      pg_clusters_TSNElogZepi_k85 == "32" ~ "MUC5AC+ goblet cell",
      pg_clusters_TSNElogZepi_k85 == "33" ~ "Basal/differentiating cell",
      pg_clusters_TSNElogZepi_k85 == "34" ~ "NK cell in epithelium",
      pg_clusters_TSNElogZepi_k85 == "35" ~ "MUC5B+ secretory cell",
      pg_clusters_TSNElogZepi_k85 == "36" ~ "FOXJ1+ ciliated cell",
      pg_clusters_TSNElogZepi_k85 == "37" ~ "CEACAM5+ peri-goblet cell"
    )
  )

metadata$CellType_Epi <- factor(metadata$CellType_Epi, levels = CellInEpi_Order)

colData(spe_epi) <- DataFrame(metadata)
metadata(spe_epi)$color_vectors$CellTypes <- CellInEpi_Color
metadata(spe_epi)$color_vectors$Histology_cell_detailed <- Histology_cell_detail_color
metadata(spe_epi)$color_vectors$Histology_cell <- Histology_cell_color
# saveRDS(spe_epi, "../input/spe_epi.rds")

# TSNE for epithelium
plot <- dittoDimPlot(spe_epi, var = "CellType_Epi",
                     reduction.use = "TSNE_log10Zepi", size = 0.05, do.label = FALSE,
                     show.axes.numbers = FALSE, ylab = "tSNE2", xlab = "tSNE1",
                     main = "", 
                     legend.title = "Cell Types in Epithelium", legend.size = 8)+ 
  scale_color_manual(values = metadata(spe_epi)$color_vectors$CellTypes, name = "Cell Types in Epithelium")  + 
  theme_void()

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
        axis.line = element_line( arrow = arrow(angle = 15, length = unit(.5, "cm"), type = "closed") )
  )

figure_layout <- c(
  area(t = 1, l = 1, b = 11, r = 11),
  area(t = 10, l = 1, b = 11, r = 2))

plot_figure <- plot + axis_plot +
  plot_layout(design = figure_layout)

# dir.create("../output/Figure_5")
pdf("../output/Figure_5/Figure_5A.pdf", height = 5, width = 7)
plot_figure
dev.off()



# ==========================================================================
# S.Figure 7A: TNSE for epithelium
# ==========================================================================
# visualize Sample ID
p1 <- dittoDimPlot(spe_epi, var = "SampleID", 
                   reduction.use = "TSNE_log10Zepi", size = 0.02, do.label = FALSE,
                   show.axes.numbers = FALSE, ylab = "tSNE2", xlab = "tSNE1",
                   main = "",
                   legend.title = "Sample ID", legend.size = 6) + 
  scale_color_manual(values = metadata(spe_epi)$color_vectors$SampleID, name = "Sample ID") + 
  theme_void() + axis_plot + plot_layout(design = figure_layout)

# dir.create("../output/SFigure_7")
pdf("../output/SFigure_7/SFigure_7A_1.pdf", height = 4, width = 5)
print(p1)
dev.off()

p2 <- dittoDimPlot(spe_epi, var = "Progression_Status", 
                   reduction.use = "TSNE_log10Zepi", size = 0.02,
                   main = "",
                   legend.title = "Progression Status", legend.size = 6) + 
  scale_color_manual(values = metadata(spe_epi)$color_vectors$Progression_Status, name = "Progression Status") +
  theme_void() + axis_plot + plot_layout(design = figure_layout)

pdf("../output/SFigure_7/SFigure_7A_2.pdf", height = 4, width = 5)
print(p2)
dev.off()

p3 <- dittoDimPlot(spe_epi, var = "Genomic_Smoking_Status", 
                   reduction.use = "TSNE_log10Zepi", size = 0.02,
                   main = "",
                   legend.title = "Smoking Status", legend.size = 6) + 
  scale_color_manual(values = metadata(spe_epi)$color_vectors$Smoking_Status, name = "Smoking Status") +
  theme_void() + axis_plot + plot_layout(design = figure_layout)

pdf("../output/SFigure_7/SFigure_7A_3.pdf", height = 4, width = 5)
print(p3)
dev.off()

p4 <- dittoDimPlot(spe_epi, var = "Histology_cell_detailed", 
                   reduction.use = "TSNE_log10Zepi", size = 0.02,
                   main = "",
                   legend.title = "Histology", legend.size = 6) + 
  scale_color_manual(values = metadata(spe_epi)$color_vectors$Histology_cell_detailed, name = "Histology") +
  theme_void() + axis_plot + plot_layout(design = figure_layout)

pdf("../output/SFigure_7/SFigure_7A_4.pdf", height = 4, width = 6)
print(p4)
dev.off()

p5 <- dittoDimPlot(spe_epi, var = "Histology_cell", 
                   reduction.use = "TSNE_log10Zepi", size = 0.02,
                   main = "",
                   legend.title = "Histology", legend.size = 6) + 
  scale_color_manual(values = metadata(spe_epi)$color_vectors$Histology_cell, name = "Broad histology") +
  theme_void() + axis_plot + plot_layout(design = figure_layout)

pdf("../output/SFigure_7/SFigure_7A_5.pdf", height = 4, width = 5)
print(p5)
dev.off()


# ==========================================================================
# Figure 5B: Heatmap of average expression of cell types in epithelium
# ==========================================================================
# compute the log-Z normalized within epithelium count matrix
count_matrix <- as.data.frame(t(assay(spe_epi, "log10")))
metadata_cell <- as.data.frame(colData(spe_epi))
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
assay(spe_epi, "log10_Z_epi") <- t(countMatrix_Log_Z1)

# saveRDS(spe_epi, "../input/spe_epi.rds")

# get the count matrix for heatmap
expMatDF <- as.data.frame(t(assay(spe_epi, "log10_Z_epi")))
expMatDF_logZstromalNormalized <- as.data.frame(t(assay(spe_epi, "log10_Z_stromalNormalized")))

# Set gene list to be displayed on heatmap
Canonical_markers_for_heatmap = c("K5",  "p40", "K8", "Ecadherin", "MUC5AC", "MUC5B",  "CEACAM5",   "SCGB1A1", "FoxJ1",
                                  "CD20", "CD3", "CD4","CD8a",  "GZMB", "PDL1", "CD68", "CD163",   "HLA_DR") 

NLRC5_target_markers_for_heatmap = c("NLRC5", "B2M",  "PSMB9", "TAP1") 

expMatDF$Cluster = spe_epi$CellType_Epi
expMatDF_logZstromalNormalized$Cluster = spe_epi$CellType_Epi

# Take average expression of feature per cluster
# canonical markers
expMatDFSummaryMean <- expMatDF %>%
  group_by(Cluster) %>%
  summarise(across(
    .cols = all_of(Canonical_markers_for_heatmap),
    .fns = list(
      mean = ~ mean(.x, na.rm = TRUE)))) %>% as.data.frame()

rownames(expMatDFSummaryMean) <- expMatDFSummaryMean$Cluster
expMatDFSummaryMean$Cluster <- NULL


# NLRC5 markers
expMatDFSummaryMean_logZstromalNormalized <- expMatDF_logZstromalNormalized %>%
  group_by(Cluster) %>%
  summarise(across(
    .cols = all_of(NLRC5_target_markers_for_heatmap),
    .fns = list(
      mean = ~ mean(.x, na.rm = TRUE)))) %>% as.data.frame()

rownames(expMatDFSummaryMean_logZstromalNormalized) <- expMatDFSummaryMean_logZstromalNormalized$Cluster
expMatDFSummaryMean_logZstromalNormalized$Cluster <- NULL

expMatDFSummaryMean = expMatDFSummaryMean[CellInEpi_Order,]
expMatDFSummaryMean_logZstromalNormalized = expMatDFSummaryMean_logZstromalNormalized[CellInEpi_Order,]

colnames(expMatDFSummaryMean) <- gsub("_mean", "", colnames(expMatDFSummaryMean))
colnames(expMatDFSummaryMean_logZstromalNormalized) <- gsub("_mean", "", colnames(expMatDFSummaryMean_logZstromalNormalized))

# cell type proportion
propTableForHeatmap <- as.matrix(table(expMatDF$Cluster)) %>% prop.table(2) *100
propTableForHeatmap = propTableForHeatmap[CellInEpi_Order,]
all.equal(rownames(expMatDFSummaryMean), names(propTableForHeatmap))

# heatmap plotting
rowAnno_right = ComplexHeatmap::rowAnnotation(
  `% Cell Types`  = ComplexHeatmap::anno_barplot(
    propTableForHeatmap, which = "row", axis = TRUE, bar_width = 0.8, width = unit(2.5, "cm"),
    gp = grid::gpar(fill = CellInEpi_Color), 
    border = TRUE,
    show_annotation_name = FALSE))


rowAnno_left = ComplexHeatmap::rowAnnotation(
  `Cell Types` = CellInEpi_Order, 
  col = list(`Cell Types` = CellInEpi_Color), 
  simple_anno_size = unit(0.3, "cm"), 
  show_legend = FALSE, 
  show_annotation_name = TRUE,
  annotation_name_rot = 90)

# Create color scale for heatmap
col1 = circlize::colorRamp2(c(-3, 0, 3),  
                            c("#1E90FF", "#FFFFFF", "#CD2626"))

col2 = circlize::colorRamp2(c(-1, 0, 2),
                            c("#02f900", "#FFFFFF", "#ffa921"))

# Create heatmap
ht1 <- ComplexHeatmap::Heatmap(matrix = expMatDFSummaryMean, 
                               heatmap_legend_param = list(
                                 col_fun = col1,
                                 title = "Average Expression\n(Log-Z normalized within epithelium)",
                                 direction = "horizontal"),
                               col = col1,
                               cluster_rows = FALSE,
                               cluster_columns = FALSE, column_names_rot = 90,
                               left_annotation = rowAnno_left, 
                               row_names_gp = grid::gpar(fontsize = 10), 
                               column_names_gp = grid::gpar(fontsize = 10), 
                               show_column_names = TRUE,
                               show_row_names = FALSE,
                               column_title = "Canonical phenotypic markers")

ht2 <- ComplexHeatmap::Heatmap(matrix = as.matrix(expMatDFSummaryMean_logZstromalNormalized),
                               heatmap_legend_param = list(
                                 col_fun = col2,
                                 title = "Average Expression\n(Log-Z stromal normalized)",
                                 direction = "horizontal"),
                               col= col2,
                               cluster_rows = FALSE,
                               cluster_columns = FALSE, column_names_rot = 90,
                               right_annotation = rowAnno_right,
                               row_names_gp = grid::gpar(fontsize = 10), 
                               column_names_gp = grid::gpar(fontsize = 10), 
                               show_column_names = TRUE,
                               show_row_names = FALSE,
                               column_title = "NLRC5 targets")

p<- ComplexHeatmap::draw(ht1+ht2, heatmap_legend_side = "bottom")

pdf("../output/Figure_5/Figure_5B.pdf", height = 4, width = 5)
print(p)
dev.off()

# ==========================================================================
# S.Figure 7B: stack plot of sample ID across cell types
# ==========================================================================
cells_df <-as.data.frame(colData(spe_epi)) %>%
  select("SampleID", "CellType_Epi", "Progression_Status", "Histology_cell") 

cells_table_df <- table(cells_df$CellType_Epi, cells_df$SampleID) %>% as.data.frame()

p <- ggplot(cells_table_df, aes(y = Freq)) +
  labs(x = "Cell Types in Epithelium", y = "Proportion (%)", title = "") + 
  geom_bar( aes(x = Var1, fill = Var2), 
            position = "fill", stat = "identity") + # 
  scale_fill_manual(values = metadata(spe_epi)$color_vectors$SampleID, name = "Sample ID") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, angle = -45, hjust = 0, vjust = 1),
    axis.title = element_text(size = 16, color = "black"),
    legend.key.height  =  unit(0.8, "lines"),
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 16),
    panel.border = element_blank(),
    panel.spacing.x = unit(1, "lines")) 

pdf("../output/SFigure_7/SFigure_7B.pdf", height = 5, width = 6)
print(p)
dev.off()



# ==========================================================================
# S.Figure 7C: stack plot of cell types across histology
# ==========================================================================

cells_table_df <- table(cells_df$CellType_Epi, cells_df$Histology_cell) %>% as.data.frame() %>%
  filter(!Var2 == "Non-annotated")

p <- ggplot(cells_table_df, aes(y = Freq)) +
  geom_bar( aes(x = Var2, fill = Var1), 
            position = "fill", stat = "identity") + # 
  labs(x = "Histology", y = "Proportion (%)", title = "") + 
  scale_fill_manual(values = metadata(spe_epi)$color_vectors$CellTypes, name = "Cell Types in Epithelium") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, angle = -45, hjust = 0, vjust = 1),
    axis.title = element_text(size = 16, color = "black"),
    legend.key.height  =  unit(0.8, "lines"),
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 16),
    panel.border = element_blank(),
    panel.spacing.x = unit(1, "lines")) 

pdf("../output/SFigure_7/SFigure_7C.pdf", height = 5, width = 6)
print(p)
dev.off()


# ==========================================================================
# Figure 5C: bubble plot for differential cell type composition across histology 
# ==========================================================================
# remove cells without histology annotation 
spe_epi_hist <- spe_epi[,!spe_epi$Histology_cell== "Non-annotated"]

metadata_epi_hist <- as.data.frame(colData(spe_epi_hist))
metadata_epi_hist$image_id <- factor(as.character(metadata_epi_hist$image_id)) 
metadata_epi_hist$Progression_Status <- factor(metadata_epi_hist$Progression_Status,levels= c("Regressive", "Progressive/ \n persistent")  )
metadata_epi_hist$Histology_cell <- factor(metadata_epi_hist$Histology_cell, levels = c("Normal epithelium", "Hyper-metaplasia", "Dysplasia"))
metadata_epi_hist <- metadata_epi_hist %>% mutate( SampleID_Hist = paste0(SampleID, "_", Histology_cell ))
metadata$SampleID_Hist <- factor(metadata_epi_hist$SampleID_Hist)
colData(spe_epi_hist) <- DataFrame(metadata_epi_hist)

# run sccomp
res_epi_fComp_histology = spe_epi_hist |>
  sccomp_estimate( 
    formula_composition = ~ Histology_cell , 
    .sample =  SampleID_Hist, 
    .cell_group = CellType_Epi, 
    bimodal_mean_variability_association = FALSE,
    mcmc_seed = 12345,
    variational_inference = FALSE,
    cores = 1  ) |>
  sccomp_remove_outliers(variational_inference = FALSE, mcmc_seed = 12345)  |>
  sccomp_test( test_composition_above_logit_fold_change = 0.1)

# saveRDS(res_epi_fComp_histology,  "../output/sccomp_CellType_Epi_fComposition_Histology.rds")

res_epi_fComp_histology1 <- res_epi_fComp_histology %>% filter(parameter != "(Intercept)") %>% 
  select(CellType_Epi, parameter, c_effect,  c_FDR, c_R_k_hat) %>%
  mutate(c_FDR_crop = ifelse(c_FDR< 0.0001, 0.0001, c_FDR))%>%
  mutate(Histology = gsub("^Histology_cell", "", parameter),
         logp = -log10(c_FDR_crop)) %>%
  mutate(sig = ifelse(is.na(c_FDR_crop) | c_FDR_crop >= 0.05, "",
                      ifelse(c_FDR_crop < 0.05 & c_FDR_crop >= 0.01, "*" ,
                             ifelse(c_FDR_crop < 0.01 & c_FDR_crop >= 0.001, "**" , "***"))) )

res_epi_fComp_histology1$Histology <- factor(res_epi_fComp_histology1$Histology, levels = c("Hyper-metaplasia", "Dysplasia") )
res_epi_fComp_histology1$CellType_Epi <- factor(res_epi_fComp_histology1$CellType_Epi, levels = rev(CellInEpi_Order))

estimate_mat <- res_epi_fComp_histology1 %>% 
  select(CellType_Epi, Histology, c_effect) %>%  
  pivot_wider(names_from = Histology, values_from = c_effect) %>% 
  as.data.frame() %>% column_to_rownames(var="CellType_Epi") %>% as.matrix()
estimate_mat = estimate_mat[CellInEpi_Order, , drop = FALSE]

logp_mat <- res_epi_fComp_histology1 %>% 
  select(CellType_Epi, Histology, logp) %>%  
  pivot_wider(names_from = Histology, values_from = logp) %>% 
  as.data.frame() %>% column_to_rownames(var="CellType_Epi") %>% as.matrix()
logp_mat = logp_mat[CellInEpi_Order, , drop = FALSE]

sig_mat <- res_epi_fComp_histology1 %>% 
  select(CellType_Epi, Histology, sig) %>%  
  pivot_wider(names_from = Histology, values_from = sig) %>% 
  as.data.frame() %>% column_to_rownames(var="CellType_Epi") %>% as.matrix()
sig_mat = sig_mat[CellInEpi_Order, , drop = FALSE]

# Create color scale for bubble plot
col_fun1 = circlize::colorRamp2(c(-1.5, 0, 1.5),
                                c("#1E90FF", "#e0e0e0", "#CD2626")) 

logp_max = max(logp_mat)
bubble_size = 3
layer_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= sqrt(pindex(logp_mat, i, j)/logp_max) * unit(bubble_size, "mm"),
              gp = gpar(fill = col_fun1(pindex(estimate_mat, i, j)), col = NA))
  grid.text(x=x, y=(y-unit(1, "mm")), pindex(sig_mat, i, j))
}

# create legend
lgd1 = ComplexHeatmap:: Legend( labels = c(1,2,3,4), row_gap = unit(1.5, "mm"), title = "-log10(FDR)",
                                graphics = list(
                                  function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(1/logp_max) * unit(bubble_size, "mm"), gp = gpar(fill = "black")),
                                  function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(2/logp_max) * unit(bubble_size, "mm"), gp = gpar(fill = "black")),
                                  function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(3/logp_max) * unit(bubble_size, "mm"), gp = gpar(fill = "black")),
                                  function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(4/logp_max) * unit(bubble_size, "mm"), gp = gpar(fill = "black")))
)

pd = packLegend(lgd1,  max_height = unit(12, "cm"), column_gap = unit(0.5, "cm"))

hp1<- ComplexHeatmap::Heatmap(matrix = estimate_mat,
                              heatmap_legend_param=list(col_fun = col_fun1, title="Estimate,\ncompared to\nnormal epithelium"), 
                              col= col_fun1,
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              left_annotation = rowAnno_left,
                              rect_gp = gpar(type = "none"),
                              layer_fun = layer_fun,
                              show_column_names = TRUE,
                              show_row_names = FALSE,
                              column_names_gp = grid::gpar(fontsize = 10), 
                              column_names_rot = 90,
                              row_names_gp = gpar(fontsize = 10),
                              border = "black"
)

p<-ComplexHeatmap::draw( hp1,  annotation_legend_list =pd) 

pdf("../output/Figure_5/Figure_5C.pdf", height = 4.5, width = 4)
print(p)
dev.off()



# ==========================================================================
# S.Figure 7D, 7E: Distribution of log-Z stromal normalized NLRC5
# ==========================================================================
# gating the epitheial cell by log-Z stromal normalized NLRC5
count_matrix <- as.data.frame(t(assay(spe_epi, "log10_Z_stromalNormalized")))
metadata <- as.data.frame(colData(spe_epi))
all.equal(rownames(metadata), rownames(count_matrix))

metadata$NLRC5_logZ_stromalNormalized = count_matrix$NLRC5
metadata$TAP1_logZ_stromalNormalized = count_matrix$TAP1
metadata$HLA_ABC_logZ_stromalNormalized = count_matrix$HLA_ABC
metadata$PSMB9_logZ_stromalNormalized = count_matrix$PSMB9
metadata$B2M_logZ_stromalNormalized = count_matrix$B2M
colData(spe_epi) <- DataFrame(metadata)


## distribution of NLRC5 expression
metadata_stats <- metadata %>% filter(!CellType_Epi %in% immuneCell_inEpi)%>%
  summarize(mean1 = mean(NLRC5_logZ_stromalNormalized) )

p<-ggplot(metadata %>% filter(!CellType_Epi %in% immuneCell_inEpi), 
          aes(x = NLRC5_logZ_stromalNormalized, y = after_stat(density), fill = SampleID)) + 
  geom_density(fill = "gray", show.legend = F) + 
  geom_density(aes(color = SampleID, fill = NULL)) +
  scale_color_manual(values = metadata(spe_epi)$color_vectors$SampleID, name = "Sample ID") +
  geom_vline(aes(xintercept = mean1), metadata_stats, color = "red", linewidth = 1) +
  labs(x = "NLRC5 (log-Z stromal normalized)", y = "Density", 
       color = "Sample ID")+
  theme_cowplot() +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key.height  =  unit(1, "lines"),
        axis.text.x = element_text(size = 14, angle = 0, vjust = 1, hjust=0.5),
        axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14)) 

print(p)

pdf("../output/SFigure_7/SFigure_7D.pdf", height = 4, width = 5.5)
print(p)
dev.off()


## NLRC5 expression on tSNE
spe_epi_rmImm <- spe_epi[, !spe_epi$CellType_Epi %in% immuneCell_inEpi]

metadata_epi_rmImm <- as.data.frame(colData(spe_epi_rmImm) )
metadata_epi_rmImm <- metadata_epi_rmImm %>% 
  mutate(NLRC5_logZ_stromalNormalized_crop = ifelse(NLRC5_logZ_stromalNormalized > 2, 2, 
                                               ifelse(NLRC5_logZ_stromalNormalized < -2, -2, 
                                                      NLRC5_logZ_stromalNormalized)))
colData(spe_epi_rmImm) <- DataFrame(metadata_epi_rmImm)

p <- dittoDimPlot(spe_epi_rmImm, var = "NLRC5_logZ_stromalNormalized_crop",  reduction.use = "TSNE_log10Zepi",
                  main = "",
                  size = 0.02, legend.show =TRUE, max=12) + 
  scale_colour_gradient2(low = "#3A3A98", high = "#832424", midpoint = 0.2, name = "Log-Z stromal\nnormalized\nexpression") + 
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) + 
  theme_void() + axis_plot + plot_layout(design = figure_layout)

pdf("../output/SFigure_7/SFigure_7E.pdf", height = 4, width = 4.5)
print(p)
dev.off()


# ==========================================================================
# Figure 5D: TSNE for epithelial cells dichotomized by NLRC5
# ==========================================================================
cutoff <- metadata_stats$mean1

metadata <- metadata %>% mutate(
  Epi_NLRC5_lowHi =  ifelse(CellType_Epi %in% immuneCell_inEpi, as.character(CellType_Epi), 
                            ifelse(NLRC5_logZ_stromalNormalized >= cutoff, "NLRC5-hi epithelial cell",
                                   ifelse(NLRC5_logZ_stromalNormalized < cutoff, "NLRC5-low epithelial cell", NA ) ) )  
)  %>% mutate(
  CellType_NLRC5_lowHi =  ifelse(CellType_Epi %in% immuneCell_inEpi,  as.character(CellType_Epi), 
                                 ifelse(NLRC5_logZ_stromalNormalized >= cutoff, paste0(CellType_Epi, "_NLRC5_hi"),
                                        ifelse(NLRC5_logZ_stromalNormalized < cutoff, paste0(CellType_Epi, "_NLRC5_low"), NA ) ) )
)

colData(spe_epi) <- DataFrame(metadata)


# set the colors vector for NLRC5 level
NLRC5Level_color <- c("#0079bc", "#99cee3")

CellTypes_NLRC5lowHi_color <- setNames(c( "#0079bc",          "#99cee3",          
                                    "#00b03a",                 "#71e9b1",           "#7D3C98",           "#D7BDE2", "gray" ),
                                 c("NLRC5-hi epithelial cell", "NLRC5-low epithelial cell", 
                                   "CD4 T cell in epithelium", "CD8 T cell in epithelium", "M1/M2 macrophage in epithelium", "NK cell in epithelium", "Undefined"))

metadata(spe_epi)$color_vectors$CellTypes_NLRC5lowHi <- CellTypes_NLRC5lowHi_color

# saveRDS(spe_epi, "../input/spe_epi.rds")

p<- dittoDimPlot(spe_epi, var = "Epi_NLRC5_lowHi",
                 reduction.use = "TSNE_log10Zepi", size = 0.05, do.label = FALSE,
                 show.axes.numbers = FALSE, ylab = "tSNE2", xlab = "tSNE1",
                 main = "",
                 legend.title = "Cell Types in Epithelium", legend.size = 8) + 
  scale_color_manual(values = metadata(spe_epi)$color_vectors$CellTypes_NLRC5lowHi, name = "Cell Types in Epithelium") +
  theme_void() + axis_plot + plot_layout(design = figure_layout)

pdf("../output/Figure_5/Figure_5D.pdf", height = 5, width = 7)
print(p)
dev.off()


# ==========================================================================
# Figure 5E: DE of TAP1, PSMB9, B2M between NLRC5 high vs low cells
# ==========================================================================
metadata_NLRC5ds <- metadata %>% filter(!Epi_NLRC5_lowHi %in% immuneCell_inEpi) %>% 
  select(Epi_NLRC5_lowHi, TAP1_logZ_stromalNormalized, PSMB9_logZ_stromalNormalized, B2M_logZ_stromalNormalized) %>% 
  pivot_longer(!Epi_NLRC5_lowHi, names_to = "Marker", values_to = "logZ_stromalNormalized") %>%
  mutate(Marker1 = gsub("_logZ_stromalNormalized", "", Marker))

metadata_NLRC5ds$Marker1 <- factor(metadata_NLRC5ds$Marker1, 
                                   levels = c("TAP1",  "PSMB9", "B2M"))

# statistics
stat.test.wilcox <- metadata_NLRC5ds %>%
  group_by(Marker1) %>%
  wilcox_test(logZ_stromalNormalized ~ Epi_NLRC5_lowHi) %>% # t_test()
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

stat.test.wilcox

p <- ggplot(metadata_NLRC5ds, aes(x = Marker1, y = logZ_stromalNormalized, fill = Epi_NLRC5_lowHi ) )+
  geom_violin( width=1, show.legend = TRUE) + 
  geom_boxplot( 
    width=0.2, color="darkgray", alpha=1, 
    position = position_dodge(width = 1),  
    outlier.color = NA, outlier.shape  = NA, show.legend = FALSE) + 
  scale_fill_manual(values = CellTypes_NLRC5lowHi_color, name = "Cell Types in Epithelium")+
  geom_signif(
    y_position = c(10, 10, 10), xmin = c(0.7, 1.7, 2.7), xmax = c(1.3, 2.3, 3.3),
    annotation = c("***", "***", "***"), tip_length = 0.02
  ) +
  labs(x = NULL, y = paste0(  "Expression level \n(log-Z stromal normalized) "),
       title = "") + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, size=14, color = "black"),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size=14),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=14),
    legend.key.height  =  unit(1, "lines"),
    legend.title = element_blank(), 
    legend.text = element_text(size=12),
    legend.position = "bottom")


pdf("../output/Figure_5/Figure_5E.pdf", height = 4, width = 5)
print(p)
dev.off()


# ==========================================================================
# S.Figure 7F: Correlation of NLRC5 targets at cell level
# ==========================================================================

# function for computing the correlation p-value
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j])
      #tmp <- cor.test(mat[, i], mat[, j], method = "spearman")
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# compute the pearson correlation and p-value
corr_data_NLRC5target = metadata %>% select(NLRC5_logZ_stromalNormalized, B2M_logZ_stromalNormalized, PSMB9_logZ_stromalNormalized, TAP1_logZ_stromalNormalized)

colnames(corr_data_NLRC5target) = c("NLRC5", "B2M", "PSMB9", "TAP1")

M_corr_NLRC5target <- cor(corr_data_NLRC5target) 
p.mat_NLRC5target <- cor.mtest(corr_data_NLRC5target)

col=rev(brewer.pal(n=10, name="RdBu")) 

pdf("../output/SFigure_7/SFigure_7F.pdf", height = 5, width = 5.5)
corrplot::corrplot(M_corr_NLRC5target, type="upper", order="original", tl.col="black", tl.srt=45, tl.cex= 1.6, is.corr = TRUE, cl.cex= 1.4, cl.ratio = 0.3,
                   p.mat = p.mat_NLRC5target, sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', pch.cex = 1.6 , pch.col = 'white', col=col, col.lim = c(-1, 1)) # order="hclust",

dev.off()


# ==========================================================================
# S.Figure 7G: differential cell type composition between NLRC5 high vs low cells
# ==========================================================================

## stack plot on left
cells_df <-as.data.frame(colData(spe_epi)) %>%
  select("SampleID", "CellType_Epi",  "Epi_NLRC5_lowHi", "Progression_Status", "Histology_cell") 

cells_table_df <- table(cells_df$CellType_Epi, cells_df$Epi_NLRC5_lowHi) %>% as.data.frame() %>%
  filter(!Var1 %in% immuneCell_inEpi, !Freq ==0) %>%
  mutate(NLRC5_level = ifelse(Var2 == "NLRC5-hi epithelial cell", "High", "Low"),
         CellTypes = factor(Var1, levels = EpiCell_Order))

p <- ggplot(cells_table_df, aes(y = Freq, x = NLRC5_level, fill = CellTypes)) +
  geom_bar( position = "fill", stat = "identity") + # 
  scale_fill_manual(values = EpiCell_Color, name = "Epitheial types") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
  labs(x = "NLRC5 levels", y = "Proportion (%)", title = "") + 
  theme_bw() + theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, vjust = 1),
    axis.title = element_text(size = 16, color = "black"),
    legend.position = "none",
    panel.border = element_blank(),
    panel.spacing.x = unit(1, "lines")) 

print(p)

pdf("../output/SFigure_7/SFigure_7G_left.pdf", height = 5, width = 2)
print(p)
dev.off()


## bubble plot in the middle and bar plot on right
spe_epi_rmImm <- spe_epi[, !spe_epi$CellType_Epi %in% immuneCell_inEpi]

metadata_epi_rmImm <- as.data.frame(colData(spe_epi_rmImm))
metadata_epi_rmImm$image_id <- factor(as.character(metadata_epi_rmImm$image_id)) 
metadata_epi_rmImm$CellType_Epi <-  factor(metadata_epi_rmImm$CellType_Epi, levels = EpiCell_Order)
metadata_epi_rmImm <- metadata_epi_rmImm %>% mutate( SampleID_LowHi = paste0(SampleID, "_", Epi_NLRC5_lowHi ))
metadata_epi_rmImm$SampleID_LowHi <- factor(metadata_epi_rmImm$SampleID_LowHi)

colData(spe_epi_rmImm) <- DataFrame(metadata_epi_rmImm)

# run sccomp to test significant
res_epi_fComp_LowHi_contrast = spe_epi_rmImm |>
  sccomp_estimate( 
    formula_composition = ~ 0 + Epi_NLRC5_lowHi, 
    .sample =  SampleID_LowHi,  
    .cell_group = CellType_Epi, 
    bimodal_mean_variability_association = FALSE,
    mcmc_seed = 12345,
    variational_inference = FALSE,
    cores = 1 
  ) |>
  sccomp_remove_outliers(variational_inference = FALSE, mcmc_seed = 12345) 

# saveRDS(res_epi_fComp_LowHi_contrast,  "../output/sccomp_CellType_Epi_fComposition_NLRC5lowHi.rds")

res_epi_fComp_LowHi_contrast <- readRDS( "../output/sccomp_CellType_Epi_fComposition_NLRC5lowHi.rds")

res_epi_fComp_LowHi_contrast_test = res_epi_fComp_LowHi_contrast |>
  sccomp_test( contrasts =  c("`Epi_NLRC5_lowHiNLRC5-hi epithelial cell` - `Epi_NLRC5_lowHiNLRC5-low epithelial cell`"), test_composition_above_logit_fold_change = 0.1)  |>
  filter(parameter != "(Intercept)") |> 
  select(CellType_Epi, parameter, c_effect,  c_FDR, c_R_k_hat) |>
  mutate(parameter = "NLRC5-hi - NLRC5-low")|>
  mutate(c_FDR_crop = ifelse(c_FDR< 0.0001, 0.0001, c_FDR),
         logp = -log10(c_FDR_crop),
         sig = ifelse(is.na(c_FDR_crop) | c_FDR_crop >= 0.05, "",
                      ifelse(c_FDR_crop < 0.05 & c_FDR_crop >= 0.01, "*" ,
                             ifelse(c_FDR_crop < 0.01 & c_FDR_crop >= 0.001, "**" , "***"))) ) |>
  mutate(CellType_Epi = factor(CellType_Epi, levels = rev(EpiCell_Order)) )


estimate_mat2 <- res_epi_fComp_LowHi_contrast_test %>% 
  select(CellType_Epi, parameter, c_effect) %>%  
  pivot_wider(names_from = parameter, values_from = c_effect) %>% 
  as.data.frame() %>% column_to_rownames(var="CellType_Epi") %>% as.matrix()
estimate_mat2 = estimate_mat2[EpiCell_Order, , drop = FALSE]

logp_mat2 <- res_epi_fComp_LowHi_contrast_test %>% 
  select(CellType_Epi, parameter, logp) %>%  
  pivot_wider(names_from = parameter, values_from = logp) %>% 
  as.data.frame() %>% column_to_rownames(var="CellType_Epi") %>% as.matrix()
logp_mat2 = logp_mat2[EpiCell_Order, , drop = FALSE]

sig_mat2 <- res_epi_fComp_LowHi_contrast_test %>% 
  select(CellType_Epi, parameter, sig) %>%  
  pivot_wider(names_from = parameter, values_from = sig) %>% 
  as.data.frame() %>% column_to_rownames(var="CellType_Epi") %>% as.matrix()
sig_mat2 = sig_mat2[EpiCell_Order, , drop = FALSE]


# Create color scale for bubble plot
col_fun2 = circlize::colorRamp2(c(-2, 0, 2),
                                c("#1E90FF", "#e0e0e0", "#CD2626")) 

logp_max = max(logp_mat2)
bubble_size = 3
layer_fun2 = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h,   gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= sqrt(pindex(logp_mat2, i, j)/4) * unit(bubble_size, "mm"),
              gp = gpar(fill = col_fun2(pindex(estimate_mat2, i, j)), col = NA))
  grid.text(x=x,y=(y-unit(1, "mm")),pindex(sig_mat2, i, j))
}

# create legend
lgd2 = ComplexHeatmap:: Legend( labels = c(1,2,3,4), row_gap = unit(1.5, "mm"), title = "-log10(FDR)",
                                graphics = list(
                                  function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(1/logp_max) * unit(bubble_size, "mm"), gp = gpar(fill = "black")),
                                  function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(2/logp_max) * unit(bubble_size, "mm"), gp = gpar(fill = "black")),
                                  function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(3/logp_max) * unit(bubble_size, "mm"), gp = gpar(fill = "black")),
                                  function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(4/logp_max) * unit(bubble_size, "mm"), gp = gpar(fill = "black")))
)


rowAnno_left = ComplexHeatmap::rowAnnotation(
  `Epithelial types` = EpiCell_Order, 
  col = list(`Epithelial types` = EpiCell_Color),
  simple_anno_size = unit(0.3, "cm"),
  show_legend = FALSE,
  show_annotation_name = TRUE,
  annotation_name_rot = 90
)

# bar plot on the right
cells_table <- table(metadata_epi_rmImm$CellType_Epi, metadata_epi_rmImm$Epi_NLRC5_lowHi) 
cells_table <- cells_table[EpiCell_Order, ]

rowAnno_right = ComplexHeatmap::rowAnnotation(
  `Number of cells`  = ComplexHeatmap::anno_barplot(cells_table, which = "row", axis = TRUE, bar_width = 0.7, width = unit(2, "cm"),
                                                    gp = grid::gpar(fill = c( "#0079bc", "#99cee3")), # 
                                                    border = FALSE,
                                                    show_annotation_name = TRUE),
  simple_anno_size = unit(0.2, "cm"),
  show_legend = FALSE)

lgd3 = ComplexHeatmap::Legend(labels = c("High", "Low"), title = "NLRC5 level", 
                              legend_gp = grid::gpar(fill = c("#0079bc", "#99cee3")))

pd = packLegend( lgd3, lgd2,  max_height = unit(11, "cm"), column_gap = unit(0.5, "cm"))

hp1<- Heatmap(estimate_mat2,
              heatmap_legend_param=list(col_fun = col_fun2, title="Estimate"), 
              col=col_fun2,
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              left_annotation = rowAnno_left,
              right_annotation = rowAnno_right,
              rect_gp = gpar(type = "none"),
              layer_fun = layer_fun2,
              show_column_names = TRUE,
              show_row_names = FALSE,
              column_names_gp = grid::gpar(fontsize = 10), 
              column_names_rot = 90,
              row_names_gp = gpar(fontsize = 10),
              border = "black"
)

p <- draw( hp1,  annotation_legend_list =pd) 

pdf("../output/SFigure_7/SFigure_7G_right.pdf", height = 5, width = 3.5)
print(p)
dev.off()


# ==========================================================================
# S.Figure 7H: stack plot of cell types between progression status
# ==========================================================================
cells_df <-as.data.frame(colData(spe_epi))
cells_table_df <- table(cells_df$SampleID, cells_df$CellType_Epi) %>% as.data.frame()
colnames(cells_table_df) <- c("Sample", "CellType",  "Freq")

m <- match(cells_table_df$Sample, cells_df$SampleID)
cells_table_df$Progression_Status <- cells_df$Progression_Status[m]


p <- ggplot(cells_table_df, aes(y = Freq)) +
  labs(x = "Sample ID", y = "Proportion (%)", title = "") + 
  facet_wrap("Progression_Status", scales = "free_x") +
  geom_bar( aes(x = Sample, fill = CellType), 
            position = "fill", stat = "identity") + # 
  scale_fill_manual(values = CellInEpi_Color, name = "Cell types in epithelium") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) +
  theme_bw() + theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 14, face = "bold"), # facet title
    strip.background = element_rect(fill = NA, color = NA),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 16, color = "black"),
    legend.key.height  =  unit(0.8, "lines"),
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 16),
    panel.border = element_blank(),
    panel.spacing.x = unit(1, "lines")) 

pdf("../output/SFigure_7/SFigure_7H.pdf", height = 5, width = 9)
print(p)
dev.off()


# ==========================================================================
# Figure 6A: Volcano plot for differential cell type composition between progression status
# ==========================================================================
res_CelltypeEpi_fComp_prog = 
  spe_epi |>
  sccomp_estimate( 
    formula_composition = ~ Progression_Status, 
    .sample = SampleID,
    .cell_group = CellType_Epi, 
    bimodal_mean_variability_association = FALSE,
    mcmc_seed = 12345,
    variational_inference = FALSE,
    cores = 1 ) |>
  sccomp_remove_outliers(variational_inference = FALSE, mcmc_seed = 12345)  |>
  sccomp_test(test_composition_above_logit_fold_change = 0.1)

# saveRDS(res_CelltypeEpi_fComp_prog,  "../output/sccomp_CellType_Epi_fComposition_Progression.rds")
res_CelltypeEpi_fComp_prog <- readRDS(  "../output/sccomp_CellType_Epi_fComposition_Progression.rds")

res_CelltypeEpi_fComp_prog1 <- res_CelltypeEpi_fComp_prog %>% 
  filter(parameter != "(Intercept)") %>% 
  select(CellType_Epi, parameter, c_effect, c_pH0, c_FDR, c_R_k_hat) 

vplot <- ggplot(res_CelltypeEpi_fComp_prog1) +
  aes(y=-log10(c_FDR), x=c_effect, 
      label =  CellType_Epi) +
  geom_point(size=5, aes(color = CellType_Epi), show.legend = FALSE) + # 
  scale_color_manual(values = CellInEpi_Color)+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", color="red", linewidth=0.5) +
  annotate("text", x=0.4, y=-log10(0.05)+0.05,
           label=paste("FDR<0.05"), size=5, fontface="bold", color = "red") +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", linewidth=0.5) +
  geom_text_repel(size = 4) +
  labs(y="-Log10(FDR)", x = "Estimate, Progressive - Regressive" ) + theme_bw()+
  theme(axis.text.x = element_text( size=16, color = "black"),
        axis.text.y = element_text(size=16, color = "black"),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18)
  ) 

# dir.create("../output/Figure_6")
pdf("../output/Figure_6/Figure_6A.pdf", height = 5, width = 6)
print(vplot)
dev.off()


# ==========================================================================
# Figure 6C: Differential NLRC5 expression between progression status across histology
# ==========================================================================
metadata1 <- colData(spe_epi) %>% as.data.frame() %>% 
  filter(!CellType_Epi %in% immuneCell_inEpi) %>% 
  filter(Histology_cell %in% c( "Normal epithelium", "Hyper-metaplasia", "Dysplasia")) 

## Mixed effect model with interaction of Progression Status and Histology
mixed_model <- lmer( NLRC5_logZ_stromalNormalized ~ Progression_Status1*Histology_cell + (1 | image_id), data = metadata1, REML=FALSE)
summary(mixed_model)


# get the emm and contrast p-value
emm_prog.his_pair = emmeans(mixed_model, specs = pairwise~ Progression_Status1 * Histology_cell, adjust = "tukey") 
emm_prog.his_contrast <- contrast(emm_prog.his_pair, interaction = "pairwise") %>% summary(infer = TRUE)
emm_prog.his_contrast

emm_df1 <- emm_prog.his_pair$emmeans %>% as.data.frame() %>% 
  mutate(Histology_cell1 = factor(Histology_cell,
                                  levels = c("Normal epithelium", "Hyper-metaplasia", "Dysplasia")))

p<-ggplot(emm_df1, aes(x=Histology_cell1, y=emmean, group=Progression_Status1, color=Progression_Status1)) +
  geom_line(aes(color=Progression_Status1), linewidth=1)+
  geom_point(aes(color=Progression_Status1), size = 3.5) +
  geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL), width=.1,  position=position_dodge(0.05)) +
  scale_color_manual(values = Progression_Status_color)+
  labs(y="Estimated Marginal Means of NLRC5", x = "Histology",
       color = "Progression Status") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, size=12, color = "black"),
        axis.text.y = element_text(size=12, color = "black"),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        text = element_text( size=12),
        legend.title = element_blank(), #element_text(size = 16), 
        legend.text  = element_text(size = 12),
        legend.position="right") 

pdf("../output/Figure_6/Figure_6C.pdf", height = 4, width = 5)
print(p)
dev.off()


# ==========================================================================
# Figure 6D: Differential NLRC5 expression between progression status across cell types and histology
# ==========================================================================
# function to compute the emm for each cell type
DiffExp_ProgHisInteractionBySubtype <- function(x) {
  metadata2 <- metadata1 %>% filter(CellType_Epi == x)  
  
  mixed_model <- lmer( NLRC5_logZ_stromalNormalized ~ Progression_Status1*Histology_cell + (1 | image_id), data = metadata2, REML=FALSE)

  emm_progHis_pair = emmeans(mixed_model, specs = pairwise~ Progression_Status1 * Histology_cell, adjust = "tukey") 
  
  emm_df <- emm_progHis_pair$emmeans %>% as.data.frame()
  
  IC_pair <- contrast(emm_progHis_pair, interaction = "pairwise") %>% summary(infer = TRUE) %>% as.data.frame() %>% mutate(CellType =as.character(x)) 
}


DiffExp_ProgHisInteractionBySubtype_list <- lapply(EpiCell_Order, DiffExp_ProgHisInteractionBySubtype)

DiffExp_ProgHisInteractionBySubtype_rbind <- rbindlist(DiffExp_ProgHisInteractionBySubtype_list, use.names=TRUE, fill=TRUE) %>% 
  select(c("CellType", "Progression_Status1_pairwise", "Histology_cell_pairwise", "estimate", "p.value"))

DiffExp_ProgHisInteractionBySubtype_rbind1 <- DiffExp_ProgHisInteractionBySubtype_rbind %>% 
  mutate(logp = -log10(p.value),
         sig = ifelse(is.na(p.value) | p.value >= 0.05, "",
                      ifelse(p.value < 0.05 & p.value >= 0.01, "*" ,
                             ifelse(p.value < 0.01 & p.value >= 0.001, "**" , "***"))) ) %>% 
  mutate(Histoloy_pair =  factor(Histology_cell_pairwise, levels = c("Normal epithelium - (Hyper-metaplasia)", "Normal epithelium - Dysplasia", "(Hyper-metaplasia) - Dysplasia"),
                                 labels = c("Hyper-metaplasia - Normal epithelium",   "Dysplasia - Normal epithelium", "Dysplasia - Hyper-metaplasia"))) %>%
  mutate( CellType = factor(CellType,  levels = rev(EpiCell_Order)))

# dot plot
cols <- c("blue", "red")

p <- ggplot(DiffExp_ProgHisInteractionBySubtype_rbind1, 
            aes(x = Histoloy_pair, y = CellType,  color = estimate,  
                label = sig, size = ifelse(is.na(logp), 1, logp)  )) + 
  geom_point(stat = "identity", na.rm = T, stroke = 2) + 
  theme_cowplot() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        strip.background = element_rect(fill = NA, color = NA),
        legend.text = element_text(size = 12), 
        axis.text.x = element_text(size = 12, angle = -30, vjust = 1, hjust=0, color = "black"),
        axis.text.y = element_text(size = 12), 
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)) + 
  labs(x = "Histology", y = "Cell type", 
       color = "EMM of NLRC5,\nProgressive - Regressive", size = "-log10(p)")  + 
  geom_text(color = "white", size = 5, vjust = 0.9) + 
  scale_color_gradient2(low = cols[1], mid = "white", high = cols[2])

print(p)


pdf("../output/Figure_6/Figure_6D.pdf", height = 5, width = 7.5)
print(p)
dev.off()

