# FDR only based on curated EC
library(tidyverse)

load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/analysis/abund0020/enzyme_trend.RData")
results_T2D <- meta_3trend

ecs <- c(
  "EC2.7.7.39: Glycerol-3-phosphate cytidylyltransferase",
  "EC2.7.8.12: CDP-glycerol glycerophosphotransferase",
  "EC2.7.8.33: N-acetylglucosaminephosphotransferase",
  
  "EC2.4.2.43: Lipid IV(A) 4-amino-4-deoxy-L-arabinosyltransferase",
  "EC2.7.7.70: D-glycero-beta-D-manno-heptose 1-phosphate adenylyltransferase",
  "EC2.4.1.129: Peptidoglycan glycosyltransferase",
  "EC2.7.1.66: Undecaprenol kinase",
  "EC2.7.1.71: Shikimate kinase",
  "EC4.2.1.10: 3-dehydroquinate dehydratase",
  "EC4.2.3.5: Chorismate synthase",
  "EC3.4.21.96: Lactocepin",
  "EC3.4.24.55: Pitrilysin",
  "EC1.2.1.12: Glyceraldehyde-3-phosphate dehydrogenase (phosphorylating)",
  "EC1.1.1.136: UDP-N-acetylglucosamine 6-dehydrogenase",
  "EC1.1.1.44: Phosphogluconate dehydrogenase (NADP(+)-dependent, decarboxylating)",
  "EC1.1.1.49: Glucose-6-phosphate dehydrogenase (NADP(+))",
  "EC3.1.6.6: Choline-sulfatase",
  "EC2.7.8.31: Undecaprenyl-phosphate glucose phosphotransferase"
)

gdata::keep(results_T2D, ecs, sure=T)

load(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/analysis/abund0020/biomarker/ecs_bio.RData")
results_mutual <- results_mutual %>% 
  filter(feature %in% ecs & k>1) %>%
  mutate(fdr=p.adjust(pval, method="BH", n=length(pval)))
results_bmi <- results_mutual[results_mutual$bio=="bmi", ]
results_ldlc <- results_mutual[results_mutual$bio=="ldlc", ]
results_hdlc <- results_mutual[results_mutual$bio=="hdlc", ]
results_tg <- results_mutual[results_mutual$bio=="tg", ]
results_homair <- results_mutual[results_mutual$bio=="homaIR", ]
results_homab <- results_mutual[results_mutual$bio=="homaB", ]
results_crp <- results_mutual[results_mutual$bio=="crp", ]

ecs_17 <- Reduce(intersect, list(results_bmi$feature,
                                 results_crp$feature,
                                 results_hdlc$feature,
                                 results_homab$feature,
                                 results_homair$feature,
                                 results_ldlc$feature,
                                 results_tg$feature))

gdata::keep(results_T2D, ecs, ecs_17,
            results_bmi, results_ldlc, results_hdlc, results_tg, results_homair,
            results_homab, results_crp, sure=T)

ecs_T2D <- results_T2D %>% filter(feature %in% ecs_17 & k>1) %>% 
  mutate(group="T2D", 
         star1=case_when(qval.fdr<0.05 ~ "#"),
         star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "*")) %>% 
  select(group, feature, coef, pval, star1, star2) %>% arrange(desc(coef)) %>% 
  mutate(cat=case_when(feature %in% c("EC2.4.2.43: Lipid IV(A) 4-amino-4-deoxy-L-arabinosyltransferase",
                                      "EC2.7.7.70: D-glycero-beta-D-manno-heptose 1-phosphate adenylyltransferase",
                                      "EC1.1.1.136: UDP-N-acetylglucosamine 6-dehydrogenase",
                                      "EC2.4.1.129: Peptidoglycan glycosyltransferase",
                                      "EC2.7.1.66: Undecaprenol kinase",
                                      "EC2.7.7.39: Glycerol-3-phosphate cytidylyltransferase",
                                      "EC2.7.8.12: CDP-glycerol glycerophosphotransferase",
                                      "EC2.7.8.33: N-acetylglucosaminephosphotransferase",
                                      "EC2.4.1.345: NO_NAME") ~ "Biosynthesis of bacterial structural components",
                       feature %in% c("EC2.7.1.71: Shikimate kinase",
                                      "EC4.2.1.10: 3-dehydroquinate dehydratase",
                                      "EC4.2.3.5: Chorismate synthase") ~ "Tryptophan biosynthesis",
                       feature %in% c("EC1.1.1.44: Phosphogluconate dehydrogenase (NADP(+)-dependent, decarboxylating)",
                                      "EC1.1.1.49: Glucose-6-phosphate dehydrogenase (NADP(+))") ~ "Glucose homeostasis",
                       feature == "EC1.2.1.12: Glyceraldehyde-3-phosphate dehydrogenase (phosphorylating)" ~ "Glucose homeostasis",
                       feature == "EC3.4.24.55: Pitrilysin" ~ "Glucose homeostasis",
                       feature == "EC3.4.21.96: Lactocepin" ~ "Catabolism of dairy components",
                       feature == "EC3.1.6.6: Choline-sulfatase"~ "Choline production",
                       feature == "EC2.7.8.31: Undecaprenyl-phosphate glucose phosphotransferase" ~ "Biosynthesis of bacterial structural components")) 
ecs_T2D$cat <- factor(ecs_T2D$cat)

# biomarker
ecs_bmi <- results_bmi %>% filter(feature %in% ecs_17) %>% 
  mutate(group="BMI", 
         star1=case_when(qval.fdr<0.05 ~ "#"),
         star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "*")) %>% 
  select(group, feature, coef, pval, star1, star2)
ecs_ldlc <- results_ldlc %>% filter(feature %in% ecs_17) %>%  
  mutate(group="LDL-C", 
         star1=case_when(qval.fdr<0.05 ~ "#"),
         star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "*")) %>% 
  select(group, feature, coef, pval, star1, star2)
ecs_hdlc <- results_hdlc %>% filter(feature %in% ecs_17) %>%  
  mutate(group="HDL-C", 
         star1=case_when(qval.fdr<0.05 ~ "#"),
         star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "*")) %>% 
  select(group, feature, coef, pval, star1, star2)
ecs_tg <- results_tg %>% filter(feature %in% ecs_17) %>%  
  mutate(group="TG", 
         star1=case_when(qval.fdr<0.05 ~ "#"),
         star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "*")) %>% 
  select(group, feature, coef, pval, star1, star2)
ecs_homair <- results_homair %>% filter(feature %in% ecs_17) %>%  
  mutate(group="HOMA-IR", 
         star1=case_when(qval.fdr<0.05 ~ "#"),
         star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "*")) %>% 
  select(group, feature, coef, pval, star1, star2)
ecs_homab <- results_homab %>% filter(feature %in% ecs_17) %>% 
  mutate(group="HOMA-B", 
         star1=case_when(qval.fdr<0.05 ~ "#"),
         star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "*")) %>% 
  select(group, feature, coef, pval, star1, star2)
ecs_crp <- results_crp %>%  filter(feature %in% ecs_17) %>% 
  mutate(group="hs-CRP", 
         star1=case_when(qval.fdr<0.05 ~ "#"),
         star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "*")) %>% 
  select(group, feature, coef, pval, star1, star2)

heatdata0 <- ecs_T2D

heatdata1 <- rbind(ecs_bmi, ecs_ldlc)
heatdata1 <- rbind(heatdata1, ecs_hdlc)
heatdata1 <- rbind(heatdata1, ecs_tg)
heatdata1 <- rbind(heatdata1, ecs_homair)
heatdata1 <- rbind(heatdata1, ecs_homab)
heatdata1 <- rbind(heatdata1, ecs_crp)
heatdata1$feature[heatdata1$feature=="EC2.4.1.345: NO_NAME"] <- "EC2.4.1.345: Phosphatidyl-myo-inositol a-mannosyltransferase"

heatdata1$group <- factor(heatdata1$group, levels=c("BMI", "HOMA-IR", "HOMA-B", "HDL-C",
                                                    "LDL-C", "TG", "hs-CRP"))
summary(heatdata0$coef)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.04545  0.11558  0.20747  0.17894  0.23000  0.44168
summary(heatdata1$coef)
# Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -0.158684 -0.026008  0.005829  0.011712  0.042856  0.291661

sum(heatdata1$coef<(-0.15)) # 1
sum(heatdata1$coef>0.3) # 1
heatdata1$coef[heatdata1$coef>=0.30] <- 0.30
heatdata1$coef[heatdata1$coef<(-0.15)] <- (-0.15)

sum(ecs_T2D$coef>0) # 14
sum(ecs_T2D$coef<0) # 3

heatdata1$feature <- gsub("EC", "EC ", heatdata1$feature)

ecs1 <- ecs_T2D[ecs_T2D$coef>0, ] %>% arrange(cat, desc(coef)) %>% pull(feature)
ecs1[ecs1=="EC2.4.1.345: NO_NAME"] <- "EC2.4.1.345: Phosphatidyl-myo-inositol a-mannosyltransferase"
ecs2 <- ecs_T2D[ecs_T2D$coef<0, ] %>% arrange(cat, coef) %>% pull(feature)

ecs1 <- gsub("EC", "EC ", ecs1)
ecs2 <- gsub("EC", "EC ", ecs2)

# cluster species
library(ComplexHeatmap)
library(circlize)

# positive
h_all <- heatdata1[,c("group", "feature", "coef")]
hall <- reshape(h_all, idvar = "feature", timevar = "group", direction = "wide")
rownames(hall) <- hall$feature
hall <- hall[,-1]
colnames(hall) <- gsub("coef.", "", colnames(hall))

col_fun = colorRamp2(c(-0.4, 0, 0.4), c("#08519c", "white", "#bd0026"))
heatall <- Heatmap(as.matrix(hall), 
                 col = col_fun,
                 cluster_columns = FALSE)
heatall <- draw(heatall)
ecs_all <- rownames(hall)[row_order(heatall)]

heatdata1$lab <- factor(heatdata1$feature, levels=rev(ecs_all))

library(ggplot2)
library(ggdendro)
library(scales)
heat1_all <- ggplot(heatdata1, aes(x=group, y=lab)) +
  geom_point(aes(fill=coef), size=11, shape=22) +
  scale_fill_gradientn(
    colors=c("#08519c","white","#bd0026"),
    values=rescale(c(-0.15,0,0.30)),
    limits=c(-0.15,0.30),
    breaks=c(-0.15,0,0.3)
  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  geom_text(aes(label=star1), color="black", size=5, nudge_y = 0.0) +
  geom_text(aes(label=star2), color="black", size=7, nudge_y = -0.2) +
  ylab("") +
  labs(fill=expression(beta * " coefficient (enzymes)")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5)) +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.8),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 60, vjust = 0, hjust = 0, color = "black", size = 12),
        axis.ticks.length=unit(0.07,"inch"),
        legend.position = "bottom") +
  guides(fill = guide_colourbar(title.position="top",
                                title.hjust = 0.5,
                                barheight = 1.2,
                                barwidth = 7))

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/figS8/fig_s8b.pdf", 
    width = 9.05, height = 7.35, onefile = F) # Open a new pdf file
print(heat1_all)
dev.off() # Close the file