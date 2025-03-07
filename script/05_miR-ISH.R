library(tidyverse) 
library(corrplot)
library(RColorBrewer)
library(ggpubr) 
library(ggsignif)
library (lme4)
library (lmerTest) 



# ==========================================================================
# Figure 4D
# ==========================================================================

### miR149 density epithelium vs non-epithelium
miR_epiVSnonEpi <- read_csv('../input/miR149_density_EpiVSnonEpi.csv')


gather_miR_epiVSnonEpi <- miR_epiVSnonEpi %>% gather("Histology", "miR149_density", miR_density_Epi, miR_density_nonEpi) %>%
  mutate(Histology1 = case_when(
    Histology == "miR_density_Epi" ~ "Epithelium",
    Histology == "miR_density_nonEpi" ~ "Non-epithelium"
  ))


p= ggplot(gather_miR_epiVSnonEpi, aes( x = Histology1, y = miR149_density)) + 
  labs(x = NULL, y = "miR-149-5p density (per area)",
       title = "" ) + 
  theme_classic() +  
  geom_boxplot(outlier.fill = NA, outlier.shape = NA, aes( fill = Histology1)) +  
  geom_line(aes(group=Sample), position = position_dodge(0.5),size=0.3) + 
  geom_point( aes( group=Sample), 
              position = position_dodge(0.5), show.legend = FALSE) + 
  ylim(0, 0.018) +
  scale_fill_manual(values = c("red", "orange")) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0, size=12),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)) +
  stat_compare_means(comparisons = list( c("Epithelium", "Non-epithelium") ), 
                     label = "p.signif", method = "wilcox.test", paired = TRUE) 


pdf("../output/Figure_4/Figure_4D.pdf", height = 4, width = 3.5)
print(p)
dev.off()


# ==========================================================================
# Figure 4E
# ==========================================================================

###  epithelial miR149 density grouped by histology in miR-ISH
miR_Histology <- read_csv("../input/miR149_density_SampleHistoLevel.csv")

miR_Histology <- miR_Histology %>% dplyr::mutate(miR_density=miR_count/nuclei_count,
              miR_density_Area=miR_count/Area)

# linear mixed model for the association of miR149 density with the interaction of progression status and histology
miR_Histology$isProgressive <- miR_Histology$Progression_Status == "Progressive/ \n persistent"
lmer_result <- lmer( miR_density ~ isProgressive*Histology_grade + (1 | Sample), input = miR_Histology, REML=FALSE)
summary(lmer_result)

# miR density adjusted for sample by by linear model
miR_Histology$miR_density.control.Sample.lm <-residuals(lm( miR_density ~ Sample, miR_Histology))

# plotting
miR_Histology$Histology_grade <- factor(miR_Histology$Histology_grade, levels = c("Normal_epithelium", "Hyperplasia_Metaplasia", "Dysplasia"),
                                            labels = c("Normal epithelium", "Hyper-metaplasia", "Dysplasia"))

progression_comparisons <- list( c("Progressive/ \n persistent", "Regressive") )
Progression_Status_color = c("#446455", "#D3DDDC")


p<- ggplot(miR_Histology, aes(x=Histology_grade, y=miR_density.control.Sample.lm,  fill=Progression_Status)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(  position = position_jitterdodge(jitter.width=0.3), size=1, alpha=1) + 
  guides(fill = guide_legend(override.aes = list(size = 0, alpha=1))) +
  geom_signif(
    y_position = c(0.5, 0.5, 0.5), xmin = c(0.7, 1.7, 2.7), xmax = c(1.3, 2.3, 3.3),
    annotation = c("ns", "ns", "*"), tip_length = 0.01
  ) +
  scale_fill_manual(values=Progression_Status_color) +
  scale_color_manual(values=Progression_Status_color) +
  labs(title="",  # miR149 density (per nucleus) adjusted for sample using linear model
       x="Histology", y = "miR-149-5p density (per nucleus) \n adjusted for sample")+
  theme_classic() + 
  theme(
    legend.position="bottom", 
    legend.title = element_blank(),
    axis.text.x = element_text(angle = -45, hjust = 0, size=12, color = "black"), #element_blank(), # 
    axis.text.y = element_text(angle = 0, hjust = 0.5, size=12, color = "black"),
    axis.title.x = element_text(size=14, color = "black"), #element_blank(),
    axis.title.y = element_text(size=12, color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())


print(p)

pdf("../output/Figure_4/Figure_4E.pdf", height = 5, width = 5)
print(p)
dev.off()







