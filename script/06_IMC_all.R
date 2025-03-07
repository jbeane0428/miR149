library(tidyverse)
library(RColorBrewer)
library(ggpubr) 
library(ggsignif)
library(ComplexHeatmap)
library (lme4)
library (lmerTest)
library(SpatialExperiment)
library(imcRtools)
library(Seurat)
library(dittoSeq)
library(patchwork)


# load spe
spe_all <- readRDS("../input/spe_all.rds")

# compute log-Z normalized count in each ROI
count_matrix <- as.data.frame(t(assay(spe_all, "log10")))
all.equal(rownames(count_matrix), colnames(spe_all))

metadata_cell <- as.data.frame(colData(spe_all))
count_matrix$ObjectID = rownames(count_matrix)
count_matrix$image_id = metadata_cell$image_id

image_ids = unique(count_matrix$image_id)

countMatrix_Log_Z = data.frame()
subset_data = c("image_id", "ObjectID")

for (id in image_ids){
  countMatrix_temp <- count_matrix %>% filter(image_id == id)
  
  celldata_subset <- countMatrix_temp[,subset_data]
  countMatrix_temp = countMatrix_temp[,!names(countMatrix_temp) %in% subset_data]
  
  countMatrix_temp <- scale(countMatrix_temp,center=T,scale=T)

  countMatrix_temp = cbind(countMatrix_temp, celldata_subset)
  
  countMatrix_Log_Z<- rbind(countMatrix_Log_Z, countMatrix_temp)
}

all.equal(countMatrix_Log_Z$ObjectID, rownames(metadata_cell))

countMatrix_Log_Z1 = countMatrix_Log_Z[,!names(countMatrix_Log_Z) %in% subset_data]

assay(spe_all, "log10_Z") <- t(countMatrix_Log_Z1)


# compute the log-Z stromal normalized count matrix
metadata_all <- as.data.frame(colData(spe_all))
count_matrix_all <- as.data.frame(t(assay(spe_all, "log10")))
all.equal(rownames(count_matrix_all), colnames(spe_all))

count_matrix_all$ObjectID = rownames(count_matrix_all)
count_matrix_all$image_id = metadata_all$image_id
count_matrix_all$isEpithelium_byCluster = metadata_all$isEpithelium_byCluster

image_ids = unique(count_matrix_all$image_id)

countMatrix_LogZ_stromalNormalized = data.frame()
subset_data = c("image_id", "ObjectID", "isEpithelium_byCluster")

for (id in image_ids){
  
  countMatrix_temp <- count_matrix_all %>% filter(image_id == id) 
  
  table_isEpi <- table(countMatrix_temp$isEpithelium_byCluster)
  N_epi = table_isEpi[["Epithelium"]]
  N_nonEpi = table_isEpi[["Non_epithelium"]]
  
  countMatrix_nonEpi <- countMatrix_temp %>% filter(isEpithelium_byCluster == "Non_epithelium")
  countMatrix_Epi <- countMatrix_temp %>% filter(isEpithelium_byCluster == "Epithelium")
  
  # sample or resample the Non-epithelial cells to make the cell proportion from epi and non-epi equal
  if (N_epi < N_nonEpi) {
    countMatrix_nonEpi1 <- countMatrix_nonEpi[sample(1:nrow(countMatrix_nonEpi), N_epi, replace=FALSE),]
  }else {
    countMatrix_nonEpi1 <- countMatrix_nonEpi[sample(1:nrow(countMatrix_nonEpi), N_epi, replace=TRUE),]
  }
  
  countMatrix_temp_stromalNormalized <- rbind(countMatrix_Epi, countMatrix_nonEpi1)
  
  celldata_subset <- countMatrix_temp[,subset_data]
  countMatrix_temp = countMatrix_temp[,!names(countMatrix_temp) %in% subset_data]
  countMatrix_temp_stromalNormalized = countMatrix_temp_stromalNormalized[,!names(countMatrix_temp_stromalNormalized) %in% subset_data]
  
  # z-score using the mean and sd from the count matrix with same epi-nonEpi proportion
  mean_stromalNormalized <- colMeans(countMatrix_temp_stromalNormalized, na.rm=TRUE)
  sd_stromalNormalized <- apply(countMatrix_temp_stromalNormalized, 2, sd)  
  countMatrix_temp_LogZ <- t((t(countMatrix_temp) - mean_stromalNormalized)/sd_stromalNormalized) 
  
  countMatrix_temp_LogZ = cbind(countMatrix_temp_LogZ, celldata_subset)
  countMatrix_LogZ_stromalNormalized<- rbind(countMatrix_LogZ_stromalNormalized, countMatrix_temp_LogZ)
}

all.equal(countMatrix_LogZ_stromalNormalized$ObjectID, rownames(metadata_all))

countMatrix_LogZ_stromalNormalized1 = countMatrix_LogZ_stromalNormalized[,!names(countMatrix_LogZ_stromalNormalized) %in% subset_data]

assay(spe_all, "log10_Z_stromalNormalized") <- t(countMatrix_LogZ_stromalNormalized1)

saveRDS(spe_all, "../input/spe_all.rds")

# ==========================================================================
# S.Figure 6A
# ==========================================================================

metadata <- as.data.frame(colData(spe_all))

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

# define color vectors for broad cell types
BroadCellType_color <- setNames(c( "#e69f00",  "#009e73" , "#56b4e9" ,  "grey"),
                                c("Epithelial Cells", "Immune Cells", "Stromal Cells", "Undefined"))

metadata(spe_all)$color_vectors$Histology_cell_detailed <- Histology_cell_detail_color
metadata(spe_all)$color_vectors$Histology_cell <- Histology_cell_color
metadata(spe_all)$color_vectors$BroadAnnotation <- BroadCellType_color


x_lab_reduc <- "UMAP1"
y_lab_reduc <- "UMAP2"

axis_plot <- ggplot(data.frame(x= 100, y = 100), aes(x = .data[["x"]], y = .data[["y"]])) +
  geom_point() +
  xlim(c(0, 10)) + ylim(c(0, 10)) +
  theme_classic() +
  ylab(y_lab_reduc) + xlab(x_lab_reduc) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(
          arrow = arrow(angle = 15, length = unit(.5, "cm"), type = "closed")
        )
  )

figure_layout <- c(
  area(t = 1, l = 1, b = 11, r = 11),
  area(t = 10, l = 1, b = 11, r = 2))


p1 <- dittoDimPlot(spe_all, var = "BroadAnnotation", 
                   reduction.use = "UMAP_log10Z_mnnForROI", size = 0.02,
                   do.label = FALSE,show.axes.numbers = FALSE,
                   main = " ",
                   legend.title = "Broad cell types", legend.size = 8) +
  scale_color_manual(values = metadata(spe_all)$color_vectors$BroadAnnotation, name = "Broad Cell Types") +
  theme_void() + axis_plot + plot_layout(design = figure_layout)


dir.create("../output/SFigure_6")
pdf("../output/SFigure_6/SFigure_6A_1.pdf", height = 4, width = 5)
print(p1)
dev.off()



p2 <- dittoDimPlot(spe_all, var = "Histology_cell_detailed", 
                   reduction.use = "UMAP_log10Z_mnnForROI", size = 0.02,
                   main = "",
                   legend.title = "Histology", legend.size = 8) + 
  scale_color_manual(values = metadata(spe_all)$color_vectors$Histology_cell_detailed, name = "Histology") +
  theme_void() + axis_plot + plot_layout(design = figure_layout)

pdf("../output/SFigure_6/SFigure_6A_2.pdf", height = 4, width = 6)
print(p2)
dev.off()


# ==========================================================================
# S.Figure 6B
# ==========================================================================


expMatDF <- as.data.frame(t(assay(spe_all, "log10_Z")))

# Set gene list to be displayed on heatmap
Canonical_markers_for_heatmap = c("K5",  "Ecadherin", "CD20", "CD3",  "HLA_DR", "aSMA", "Vimentin") 

# Add cluster information to data frame
expMatDF$Cluster = spe_all$BroadAnnotation

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

#Change feature names to names that are may be more familiar with audience
colnames(expMatDFSummaryMean) <- gsub("_mean", "", colnames(expMatDFSummaryMean))

# bar plot for cell type proportion
propTableForHeatmap4 <- as.matrix(table(expMatDF$Cluster)) %>% prop.table(2) *100
propTableForHeatmap4 <- propTableForHeatmap4[rownames(expMatDFSummaryMean),]

rowAnno_right = ComplexHeatmap::rowAnnotation(
  `% Cell Types`  = ComplexHeatmap::anno_barplot(propTableForHeatmap4, which = "row", axis = TRUE, bar_width = 0.8, width = unit(2.5, "cm"),
                                                 gp = grid::gpar(fill = c( "#e69f00",  "#009e73" , "#56b4e9" ,  "grey")), # 
                                                 border = TRUE,
                                                 show_annotation_name = FALSE))


rowAnno_left = ComplexHeatmap::rowAnnotation(
  `Cell Types` = c("Epithelial Cells", "Immune Cells", "Stromal Cells", "Undefined"), 
  col = list(`Cell Types` = 
               setNames(c( "#e69f00",  "#009e73" , "#56b4e9" ,  "grey"),
                        c("Epithelial Cells", "Immune Cells", "Stromal Cells", "Undefined"))), 
  simple_anno_size = unit(0.3, "cm"),
  show_legend = FALSE
)


# Create color scale for heatmap
col1 = circlize::colorRamp2(c(-3, 0, 3),
                            c("#1E90FF", "#FFFFFF", "#CD2626"))

# Create heatmap
ht1 <- ComplexHeatmap::Heatmap(matrix = expMatDFSummaryMean, 
                               heatmap_legend_param = list(
                                 col_fun = col1,
                                 title = "Average Expression\n(Log-Z normalized)",
                                 direction = "horizontal"),
                               col = col1,
                               cluster_rows = FALSE,
                               cluster_columns = FALSE, column_names_rot = 90,
                               left_annotation = rowAnno_left, 
                               right_annotation = rowAnno_right,
                               row_names_gp = grid::gpar(fontsize = 10), 
                               column_names_gp = grid::gpar(fontsize = 10), 
                               show_column_names = TRUE,
                               show_row_names = FALSE,
                               column_title = "Canonical markers")


p <- ComplexHeatmap::draw(ht1, heatmap_legend_side = "bottom") 

pdf("../output/SFigure_6/SFigure_6B.pdf", height = 4, width = 5)
print(p)
dev.off()




