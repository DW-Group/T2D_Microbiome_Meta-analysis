setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/figS7")
library(tidyverse)

# newly T2D in SOL
#######################################
#     check incident T2D in SOL
#######################################
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/1.preprocess/SOL/species_sol_excl.RData")
newT2D_SOL <- metadata2881 %>% filter(DIABETES5 !=3 & DIABETES5_V2==3) %>% pull(burk_id)
gdata::keep(newT2D_SOL, sure=TRUE)

load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/analysis/abund0020/species_T2D_main.RData")
gdata::keep(newT2D_SOL, meta_3, meta3, meta_1, meta1, 
            species_T2D_zero, metadata_T2D, sure=TRUE)

metadata_newT2D <- metadata_T2D %>% 
  filter(study=="Wu_2020a (SWE)" |
           study=="Wu_2020b (SWE)" |
           study=="Zhong_2019 (CHN)" |
           study=="SOL (USA)" & (status_new=="Con" | id %in% newT2D_SOL)) %>% 
  group_by(study) %>% filter(!is.na(bmi)) %>% 
  mutate(bmiq=cut(bmi, breaks=c(quantile(bmi, probs = seq(0, 1, by = 0.2))), 
                  labels = F, include.lowest = TRUE))

library(fastDummies)
metadata_newT2D <- dummy_cols(metadata_newT2D, 
                              select_columns = "bmiq")
metadata_newT2D <- as.data.frame(metadata_newT2D)
rownames(metadata_newT2D) <- metadata_newT2D$id

metadata_newT2D$study <- as.factor(as.character(metadata_newT2D$study))
species_newT2D_zero <- species_T2D_zero[, metadata_newT2D$id]

library(MMUPHin)
meta33 <- lm_meta(feature_abd = as.matrix(species_newT2D_zero),
                  batch = "study",
                  exposure = "status_new",
                  covariates = c("age", "sex", "bmiq_2", "bmiq_3", "bmiq_4", "bmiq_5"),
                  data = metadata_newT2D,
                  control = list(verbose = T,
                                 rma_method = "FE",
                                 transform = "LOG",
                                 output = "species_newT2D"))
meta_33 <- meta33$meta_fits %>% filter(!is.na(coef) & k>1)
sum(meta_33$qval.fdr<0.25) # 8
sum(meta_33$qval.fdr<0.1) # 0

meta_1v2 <- meta_1 %>% filter(feature %in% meta_33$feature)
meta_3v2 <- meta_3 %>% filter(feature %in% meta_33$feature)
identical(meta_3v2$feature, meta_33$feature) # TRUE

cor(meta_1v2$coef, meta_33$coef, method="spearman") # 0.68
cor(meta_3v2$coef, meta_33$coef, method="spearman") # 0.81
plot(meta_3v2$coef, meta_33$coef)

#######################################
#        compare among SOL
#######################################
allT2D <- meta3$maaslin_fits$`SOL (USA)` %>% filter(!is.na(coef))
newT2D <- meta33$maaslin_fits$`SOL (USA)` %>% filter(!is.na(coef))
identical(allT2D$feature, newT2D$feature)
cor(allT2D$coef, newT2D$coef) # 0.77
cor(allT2D$coef, newT2D$coef, method="spearman") # 0.76

compareSOL <- data.frame(feature=allT2D$feature,
                         allT2D=allT2D$coef,
                         newT2D=newT2D$coef)

library(ggplot2)
p1 <- ggplot(aes(x=allT2D, y=newT2D), data=compareSOL) +
  geom_point() +
  geom_hline(yintercept = 0, col="red") +
  geom_vline(xintercept = 0, col="red") +
  labs(x="Coefficients for all T2D",
       y="Coefficients for\nnewly-diagnosed T2D",
       title="") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(colour = "black")) +
  annotate(geom="text", x=-0.5, y=0.88, col="black", label="r=0.76")

gdata::keep(compareSOL, p1, sure=T)

load(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/analysis/abund0020/metf_species/species_nometf_noinsul.RData")
identical(meta_3v2$feature, meta_noinsul_results$feature)
cor(meta_3v2$coef, meta_noinsul_results$coef, method="spearman") # 0.90

compareINS <- data.frame(x=meta_3v2$coef,
                         y=meta_noinsul_results$coef)

p2 <- ggplot(aes(x=x, y=y), data=compareINS) +
  geom_point() +
  geom_hline(yintercept = 0, col="red") +
  geom_vline(xintercept = 0, col="red") +
  labs(x="Coefficients for all T2D",
       y="Coefficients for\nT2D without insulin use",
       title="") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(colour = "black")) +
  annotate(geom="text", x=-0.65, y=0.72, col="black", label="r=0.90")

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/figS7/fig_s7.pdf",
    width = 6, height = 3, onefile = F) # Open a new pdf file
egg::ggarrange(p1, p2, nrow=1, ncol=2)
dev.off() # Close the file