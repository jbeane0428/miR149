library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggpubr) 
library(ggsignif)
library(rstatix)
library(ComplexHeatmap)
library (lme4)
library (lmerTest)
library(SpatialExperiment)
library(imcRtools)
library(Seurat)
library(dittoSeq)
library(patchwork)
library(BiocParallel)


# Progression status
Regressive = c("001-0482", "1-1117", "11-433", "25-383", "29-379", "29-624", "29-625", "11-340", "25-541")
Progressive = c( "14-394", "8-838", "21-319", "21-582", "22-878",  "15-570", "1-481", "21-720", "15-402", "21-1089")

# combine the immune cell types annotation for epithelium and non-epithelium into spe_all metadata
spe_all <- readRDS( "../input/spe_all.rds")
spe_epi <- readRDS( "../input/spe_epi.rds")
spe_stroma <- readRDS( "../input/spe_stroma.rds")

metadata_all <- as.data.frame(colData(spe_all)) 

cellAnnotation_epi <- as.data.frame(colData(spe_epi)) %>% 
  select(uniqueCellID,   CellType_Epi, Epi_NLRC5_lowHi, CellType_NLRC5_lowHi) %>%
  dplyr::rename(CellType = CellType_Epi)

cellAnnotation_stroma <- as.data.frame(colData(spe_stroma)) %>% 
  select(uniqueCellID,   CellType_Stroma) %>% 
  mutate(Epi_NLRC5_lowHi = CellType_Stroma,
         CellType_NLRC5_lowHi = CellType_Stroma) %>%
  dplyr::rename(CellType = CellType_Stroma)

cellAnnotation_combined = rbind(cellAnnotation_epi, cellAnnotation_stroma)

metadata_all <- merge(metadata_all, cellAnnotation_combined, by = c("uniqueCellID"), all.x = TRUE)

metadata_all <- metadata_all %>% mutate(
  Epi_NLRC5_lowHi =  ifelse(is.na(Epi_NLRC5_lowHi), "Undefined", as.character(Epi_NLRC5_lowHi)),
  Epi_NLRC5_lowHi_combined = case_when(
    Epi_NLRC5_lowHi == "NLRC5-low epithelial cell"   ~   "NLRC5-low epithelial cell",                      
    Epi_NLRC5_lowHi == "NLRC5-hi epithelial cell" ~ "NLRC5-hi epithelial cell",
    Epi_NLRC5_lowHi == "NK cells/Neutrophil" ~ "NK cell/Neutrophil",
    Epi_NLRC5_lowHi == "NK cell in epithelium" ~ "NK cell/Neutrophil",
    Epi_NLRC5_lowHi == "CD8 T cells"  ~   "CD8 T cell",
    Epi_NLRC5_lowHi == "CD4 T cells"  ~  "CD4 T cell",  
    Epi_NLRC5_lowHi == "CD8 T cell in epithelium"  ~   "CD8 T cell", 
    Epi_NLRC5_lowHi == "CD4 T cell in epithelium"  ~   "CD4 T cell",
    Epi_NLRC5_lowHi == "CD8 T/CD4 T/M1 macrophage cells" ~ "CD8 T/CD4 T/M1 macrophage",
    Epi_NLRC5_lowHi == "B cells" ~ "B cells",
    Epi_NLRC5_lowHi == "M1/M2 macrophages" ~  "M1/M2 macrophage",
    Epi_NLRC5_lowHi == "M1/M2 macrophage in epithelium" ~  "M1/M2 macrophage",  
    Epi_NLRC5_lowHi == "M1 macrophages" ~  "M1 macrophage", 
    Epi_NLRC5_lowHi == "M2 macrophages" ~  "M2 macrophage",
    Epi_NLRC5_lowHi == "Vim+ stroma cells" ~  "VIM+ stromal cell", 
    Epi_NLRC5_lowHi == "SMA+ stroma cells" ~ "SMA+ stromal cell",
    Epi_NLRC5_lowHi == "SMA+Vim+ stroma cells" ~  "SMA+VIM+ stromal cell",
    Epi_NLRC5_lowHi == "Undefined stroma"  ~   "Undefined",
    Epi_NLRC5_lowHi == "Undefined"  ~   "Undefined",
    Epi_NLRC5_lowHi == ""  ~   "Undefined"  ), 
  CellType_NLRC5_lowHi_combined = ifelse(Epi_NLRC5_lowHi_combined %in% c("NLRC5-low epithelial cell", "NLRC5-hi epithelial cell"),
                                         CellType_NLRC5_lowHi, Epi_NLRC5_lowHi_combined),
  CellType_combined = ifelse(Epi_NLRC5_lowHi_combined %in% c("NLRC5-low epithelial cell", "NLRC5-hi epithelial cell"),
                             as.character(CellType), Epi_NLRC5_lowHi_combined) )

rownames(metadata_all) <- metadata_all$uniqueCellID
colData(spe_all) <- DataFrame(metadata_all)

# ==========================================================================
# Figure 5I: spatial proximity between NLRC5 high/low cells and immune cells
# ==========================================================================
spe_all <- buildSpatialGraph(spe_all, img_id = "image_id", type = "expansion", threshold = 30) 

# test the spatial proximity based on NLRC5 hi/low epithelial cells annotation and combined immune annotation
out_ImmuneMerge <- testInteractions(spe_all, 
                        group_by = "SampleID",
                        label = "Epi_NLRC5_lowHi_combined", 
                        colPairName = "expansion_interaction_graph", 
                        iter = 500,
                        BPPARAM = SerialParam(RNGseed = 221029))

saveRDS(out_ImmuneMerge, "../output/spatial_NLRC5LowHi_centroidExpansion30.rds")


# compute the average of spatial proximity difference between NLRC5-high and NLRC5-low epithelial cells with respect to various immune types 
proximity_compare <- out_ImmuneMerge %>% as_tibble() %>% 
  dplyr::mutate(
    SampleID = as.vector(str_extract_all(group_by, "[\\d]+\\-[\\d]+", simplify = TRUE)) )%>%
  select( from_label, to_label, SampleID, p_lt) %>%
  filter( grepl('CD4|CD8|macrophage|NK|B', from_label) ) %>% 
  filter( grepl('NLRC5', to_label) ) %>%
  mutate(NLRC5_level = str_extract(to_label, "NLRC5-hi|NLRC5-low")) %>%
  dplyr::rename(ImmuneCellType = from_label)

proximity_compare_avg <- proximity_compare %>% 
  group_by(ImmuneCellType, NLRC5_level) %>%
  summarize(avg_proximity_score = mean(p_lt, na.rm = TRUE)) %>%
  pivot_wider(names_from = NLRC5_level, values_from = avg_proximity_score) %>%
  mutate(proximity_change = `NLRC5-hi` - `NLRC5-low`) 
  
# compute the p-value  
proximity_compare_wilcox <- proximity_compare %>% group_by(ImmuneCellType) %>%
  wilcox_test(p_lt ~ NLRC5_level, paired = TRUE) %>% 
  select(ImmuneCellType, p) 


proximity_compare_avg_wilcox <-  merge(proximity_compare_avg, proximity_compare_wilcox,  
                                      by = c("ImmuneCellType") ,  all.x = TRUE)  %>%
  mutate(ImmuneCellType = forcats::fct_reorder(.f = ImmuneCellType,
                                               .x = proximity_change,
                                               .desc = TRUE)) %>% 
  mutate(`Significant change` = ifelse(p <0.055, "Yes", "No"))


p1 <-  ggplot(proximity_compare_avg_wilcox, aes(x = ImmuneCellType, y = proximity_change, fill= `Significant change`)) +
  geom_bar(stat="identity") + 
  labs(x = "Immune Cell Types", y = "Changes of the spatial proximity of\n immune cell types and epithelial cells \n(NLRC5 hi - NLRC5 low)") + 
  theme_classic() + theme(
    axis.text = element_text( size = 12, color = "black"),
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    legend.title = element_text(size = 12), 
    legend.text  = element_text(size = 12),
    legend.position="bottom") 

p1

pdf("../output/Figure_5/Figure_5I.pdf", height = 6, width = 4.5)
print(p1)
dev.off()


# ==========================================================================
# Figure 5J: spatial proximity between NLRC5 high/low epithelial subtypes and CD8 T cells
# ==========================================================================
out_CellType_NLRC5LowHi <- testInteractions(spe_all, 
                                                     group_by = "SampleID",
                                                     label = "CellType_NLRC5_lowHi_combined", 
                                                     colPairName = "expansion_interaction_graph", # "neighborhood",
                                                     iter = 500,
                                                     BPPARAM = SerialParam(RNGseed = 221029))

# saveRDS(out_CellType_NLRC5LowHi, "../output/spatial_CellType_NLRC5LowHi_centroidExpansion30.rds")

proximity_compare2 <- out_CellType_NLRC5LowHi %>% as_tibble() %>% 
  dplyr::mutate(
    SampleID = as.vector(str_extract_all(group_by, "[\\d]+\\-[\\d]+", simplify = TRUE)) )%>%
  select( from_label, to_label, SampleID, p_lt) %>%
  filter( grepl('CD8 T cell', from_label) ) %>% 
  filter( grepl('NLRC5', to_label) ) %>%
  mutate(NLRC5_level = str_extract(to_label, "NLRC5_hi|NLRC5_low"),
         EpiCellType = gsub("_NLRC5_hi$|_NLRC5_low$", "", to_label))

proximity_compare2_avg <- proximity_compare2 %>% 
  group_by(from_label, to_label) %>%
  summarize(avg_proximity_score = mean(p_lt, na.rm = TRUE)) %>%
  mutate(NLRC5_level = str_extract(to_label, "NLRC5_hi|NLRC5_low"),
         EpiCellType = gsub("_NLRC5_hi$|_NLRC5_low$", "", to_label)) %>%
  select(-to_label) %>%
  pivot_wider(names_from = NLRC5_level, values_from = avg_proximity_score) %>%
  mutate(proximity_change = NLRC5_hi - NLRC5_low) 


proximity_compare2_wilcox <- proximity_compare2 %>% group_by(EpiCellType) %>%
  wilcox_test(p_lt ~ NLRC5_level, paired = TRUE) %>% 
  select(EpiCellType, p) 


proximity_compare2_avg_wilcox <-  merge(proximity_compare2_avg, proximity_compare2_wilcox,  # joined_df3
                                      by = c("EpiCellType") ,  all.x = TRUE, all.y = FALSE)  %>%
  mutate(EpiCellType = forcats::fct_reorder(.f = EpiCellType,
                                            .x = proximity_change,
                                            .desc = TRUE)) %>% 
  mutate(`Significant change` = ifelse(p <0.05, "Yes", "No"))


p2 <-  ggplot(proximity_compare2_avg_wilcox, aes(x = EpiCellType, y = proximity_change, fill= `Significant change`)) +
  geom_bar(stat="identity") + 
  labs(x = "Epithelial Cell Subtypes", y = "Changes of the spatial proximity of\n CD8 T and epithelial subtypes\n(NLRC5 hi - NLRC5 low)") + 
  theme_classic() + theme(
    axis.text = element_text( size = 12, color = "black"),
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 1),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    legend.title = element_text(size = 12), 
    legend.text  = element_text(size = 12),
    legend.position="bottom") 

p2

pdf("../output/Figure_5/Figure_5J.pdf", height = 6.5, width = 5.5)
print(p2)
dev.off()


# ==========================================================================
# S. Figure 8D: change of spatial proximity between NLRC5 high/low epithelial subtypes and immune cells
# ==========================================================================
ImmCells_combined <- c( "NK cell/Neutrophil", "CD8 T cell",  "M1 macrophage", "M1/M2 macrophage", 
               "CD4 T cell", "M2 macrophage", "CD8 T/CD4 T/M1 macrophage",  "B cells")

proximity_compare3 <- out_CellType_NLRC5LowHi %>% as_tibble() %>% 
  dplyr::mutate(
    SampleID = as.vector(str_extract_all(group_by, "[\\d]+\\-[\\d]+", simplify = TRUE)) )%>%
  select( from_label, to_label, SampleID, p_lt) %>%
  filter( from_label %in% ImmCells_combined  ) %>% 
  filter( grepl('NLRC5', to_label) ) %>%
  mutate(NLRC5_level = str_extract(to_label, "NLRC5_hi|NLRC5_low"),
         EpiCellType = gsub("_NLRC5_hi$|_NLRC5_low$", "", to_label)) %>%
  dplyr::rename(ImmuneCellType = from_label)

proximity_compare3_avg <- proximity_compare3 %>% 
  group_by(ImmuneCellType, EpiCellType, NLRC5_level) %>%
  summarize(avg_proximity_score = mean(p_lt, na.rm = TRUE)) %>%
  pivot_wider(names_from = NLRC5_level, values_from = avg_proximity_score) %>%
  mutate(proximity_change = NLRC5_hi - NLRC5_low) 


proximity_compare3_wilcox <- proximity_compare3 %>% group_by(ImmuneCellType, EpiCellType) %>%
  wilcox_test(p_lt ~ NLRC5_level, paired = TRUE) %>% 
  select(ImmuneCellType, EpiCellType, p) 

  
proximity_compare3_avg_wilcox <-  merge(proximity_compare3_avg, proximity_compare3_wilcox,  
                                        by = c("ImmuneCellType", "EpiCellType") ,  all.x = TRUE, all.y = FALSE)  %>%
    mutate(logp = -log10(p),
           sig = ifelse(is.na(p) | p >= 0.05, "",
                        ifelse(p < 0.05 & p >= 0.01, "*" ,
                               ifelse(p < 0.01 & p >= 0.001, "**" , "***"))) )
  
proximity_compare3_avg_wilcox$ImmuneCellType <- factor(proximity_compare3_avg_wilcox$ImmuneCellType, levels =
                                                      c("CD8 T cell", "M2 macrophage",  "M1 macrophage", "B cells", "M1/M2 macrophage", "NK cell/Neutrophil", "CD4 T cell", "CD8 T/CD4 T/M1 macrophage"))
  
proximity_compare3_avg_wilcox$EpiCellType <- factor(proximity_compare3_avg_wilcox$EpiCellType, levels =
                                                      rev(c(   "Basal cell",              
                                                               "Basal/differentiating cell",
                                                               "KRT5 low basal cell",
                                                               "MUC5AC+ goblet cell",
                                                               "MUC5B+ secretory cell",
                                                               "MUC5B low secretory cell",
                                                               "CEACAM5+ peri-goblet cell",                      
                                                               "CEACAM5+SCGB1A1+ secretory cell",                                    
                                                               "FOXJ1+ ciliated cell")))
  
  cols <- c("blue", "red")
  p <- ggplot(proximity_compare3_avg_wilcox, aes(x = ImmuneCellType, y = EpiCellType,  color = proximity_change,  label = sig,
                                               size = ifelse(is.na(logp), 1, logp), shape = ifelse(is.na(proximity_change), "Missing", "Present"))) + 
    geom_point(stat = "identity", na.rm = T, stroke = 2) + 
    theme_cowplot() +
    theme(axis.text = element_text(size = 12), 
          legend.text = element_text(size = 10), 
          axis.text.x = element_text(angle = -60, vjust = 1, hjust=0),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          legend.position = "right", 
          legend.box="vertical") + 
    labs(x = "Immune Cell Types", y = "Epithelial Cell Subtypes", 
         title = "",
         color = "Change in spatial proximity\n(NLRC5 hi - NLRC5 low)", size = "-log10(p)") + guides( shape = F) +  # color = F
    geom_text(color = "white", size = 5, vjust = 0.9) + 
    scale_color_gradient2(low = cols[1], mid = "white", high = cols[2])
  
  print(p)
  

pdf("../output/SFigure_8/SFigure_8D.pdf", height = 6, width = 9)
print(p)
dev.off()
