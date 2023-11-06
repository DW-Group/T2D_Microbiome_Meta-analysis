library(tidyverse)

# only biomarker, no T2D
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/figure3a.RData")
gdata::keep(heatmap1, heatmap2,sure=T)

load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig3/figure3a_spe.RData")

load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/biomarker/spe_bio_2SD.RData")
results_mutual <- results_mutual %>% 
  filter(feature %in% spe_all & k>1) %>%
  mutate(fdr=p.adjust(pval, method="BH", n=length(pval)))

results_bmi <- results_mutual[results_mutual$bio=="bmi", ]
results_ldlc <- results_mutual[results_mutual$bio=="ldlc", ]
results_hdlc <- results_mutual[results_mutual$bio=="hdlc", ]
results_tg <- results_mutual[results_mutual$bio=="tg", ]
results_homair <- results_mutual[results_mutual$bio=="homaIR", ]
results_homab <- results_mutual[results_mutual$bio=="homaB", ]
results_crp <- results_mutual[results_mutual$bio=="crp", ]

spe_42 <- Reduce(intersect, list(results_bmi$feature,
                                 results_crp$feature,
                                 results_hdlc$feature,
                                 results_homab$feature,
                                 results_homair$feature,
                                 results_ldlc$feature,
                                 results_tg$feature))
spe_1 <- intersect(heatmap1$feature, spe_42)
spe_2 <- intersect(heatmap2$feature, spe_42)

gdata::keep(results_bmi, results_ldlc, results_hdlc, results_tg, results_homair,
            results_homab, results_crp,
            spe_42, spe_1, spe_2, sure=T)

# biomarker
spe_bmi <- results_bmi %>% filter(feature %in% spe_42) %>% mutate(group="BMI", 
                                                                  star=case_when(fdr<0.25 ~ "*",
                                                                                 TRUE ~ "")) %>% 
  select(group, feature, coef, pval, star)
spe_ldlc <- results_ldlc %>% filter(feature %in% spe_42) %>% mutate(group="LDL-C", 
                                                                    star=case_when(fdr<0.25 ~ "*",
                                                                                   TRUE ~ "")) %>% 
  select(group, feature, coef, pval, star)
spe_hdlc <- results_hdlc %>% filter(feature %in% spe_42) %>% 
  mutate(group="HDL-C", 
         star=case_when(fdr<0.25 ~ "*",
                        TRUE ~ "")) %>% 
  select(group, feature, coef, pval, star)
spe_tg <- results_tg %>% filter(feature %in% spe_42) %>% 
  mutate(group="TG", 
         star=case_when(fdr<0.25 ~ "*",
                        TRUE ~ "")) %>% 
  select(group, feature, coef, pval, star)
spe_homair <- results_homair %>% filter(feature %in% spe_42) %>% 
  mutate(group="HOMA-IR", 
         star=case_when(fdr<0.25 ~ "*",
                        TRUE ~ "")) %>% 
  select(group, feature, coef, pval, star)
spe_homab <- results_homab %>% filter(feature %in% spe_42) %>% 
  mutate(group="HOMA-B", 
         star=case_when(fdr<0.25 ~ "*",
                        TRUE ~ "")) %>% 
  select(group, feature, coef, pval, star)
spe_crp <- results_crp %>% filter(feature %in% spe_42) %>% 
  mutate(group="hs-CRP", 
         star=case_when(fdr<0.25 ~ "*",
                        TRUE ~ "")) %>% 
  select(group, feature, coef, pval, star)

heatdata1 <- rbind(spe_bmi, spe_ldlc)
heatdata1 <- rbind(heatdata1, spe_hdlc)
heatdata1 <- rbind(heatdata1, spe_tg)
heatdata1 <- rbind(heatdata1, spe_homair)
heatdata1 <- rbind(heatdata1, spe_homab)
heatdata1 <- rbind(heatdata1, spe_crp)

heatdata1$lab <- gsub("s__", "", heatdata1$feature)
heatdata1$lab <- gsub("_", " ", heatdata1$lab)

heatdata1$group <- factor(heatdata1$group, levels=c("BMI", "HOMA-IR", "HOMA-B", "HDL-C",
                                                    "LDL-C", "TG", "hs-CRP"))
summary(heatdata1$coef)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.66397 -0.13523 -0.01094 -0.01611  0.09241  0.76419

# cluster species
library(ComplexHeatmap)
library(circlize)
library(ggdendro)

# all
h_all <- heatdata1[,c("group", "lab", "coef")]
hall <- reshape(h_all, idvar = "lab", timevar = "group", direction = "wide")
rownames(hall) <- hall$lab
hall <- hall[,-1]
colnames(hall) <- gsub("coef.", "", colnames(hall))

col_fun = colorRamp2(c(-0.8, 0, 0.8), c("#08519c", "white", "#bd0026"))
heatall <- Heatmap(as.matrix(hall), 
                 col = col_fun,
                 cluster_columns = FALSE)
heatall <- draw(heatall)
spe_all <- rownames(hall)[row_order(heatall)]
spe_all_dend <- row_dend(heatall)
spe_all_dend_data <- dendro_data(spe_all_dend, type = "rectangle")

heatdata1$lab <- factor(heatdata1$lab, levels=rev(spe_all))

library(ggplot2)
library(scales)
heat_all <- ggplot(segment(spe_all_dend_data)) + 
  geom_segment(aes(x = -y, y = -x, xend = -yend, yend = -xend), size=0.8) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

heat1_all <- ggplot(heatdata1, aes(x=group, y=lab)) +
  geom_point(aes(fill=coef), size=11, shape=22) +
  scale_fill_gradient2(low="#08519c",
                       mid="white",
                       high="#bd0026",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.8,0,0.8),
                       limits=c(-0.8,0.8)) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  geom_text(aes(label=star), color="black", size=7, nudge_y = -0.2) +
  labs(color=expression(beta *" coefficients")) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.8),
        axis.text.y = element_text(size = 12, face = "italic", colour = "black"),
        axis.text.x = element_text(angle = 60, vjust = 0, hjust = 0, color = "black", size = 12),
        axis.ticks.length=unit(0.07,"inch"),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 1.5))

blank <- ggplot()+theme_void()

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/biomarker/figure_species_v4.pdf", 
    width = 5.6, height = 14.2, onefile = F) # Open a new pdf file
egg::ggarrange(heat_all, blank, heat1_all, 
               nrow=1, ncol=3, widths=c(2,-0.7,7))
dev.off() # Close the file

save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/biomarker/figure_species_v4.RData")
