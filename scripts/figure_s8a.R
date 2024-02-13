library(tidyverse)

load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/pathway_T2D.RData")
results_T2D <- meta_3

pwy <- c(
  "GLYCOLYSIS",
  "HEXITOLDEGSUPER-PWY",
  "METHGLYUT-PWY",
  "PANTO-PWY",
  "PENTOSE-P-PWY",
  "SER-GLYSYN-PWY",
  "FASYN-ELONG-PWY",
  "PWY-5347",
  "PWY-5121",
  "PWY-5971",
  "PWY-6901",
  "PWY-6168",
  "PWY-7357",
  "SULFATE-CYS-PWY",
  "TEICHOICACID-PWY",
  "NAGLIPASYN-PWY",
  "LACTOSECAT-PWY"
)

pwy <- intersect(results_T2D$feature, pwy)

gdata::keep(results_T2D, pwy, sure=T)

load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/biomarker/pwy_bio_2SD.RData")

results <- results %>% 
  filter(feature %in% pwy & k>1) %>%
  mutate(fdr=p.adjust(pval, method="BH", n=length(pval)))

results_bmi <- results[results$bio=="bmi", ]
results_ldlc <- results[results$bio=="ldlc", ]
results_hdlc <- results[results$bio=="hdlc", ]
results_tg <- results[results$bio=="tg", ]
results_homair <- results[results$bio=="homaIR", ]
results_homab <- results[results$bio=="homaB", ]
results_crp <- results[results$bio=="crp", ]

gdata::keep(results_T2D, pwy,
            results_bmi, results_ldlc, results_hdlc, results_tg, results_homair,
            results_homab, results_crp, sure=T)

pwy_T2D <- results_T2D %>% filter(feature %in% pwy) %>% 
  mutate(group="T2D", 
         star1=case_when(qval.fdr<0.05 ~ "+"),
         star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "#"),
         star3=case_when(qval.fdr>=0.1 &  qval.fdr<0.25 ~ "*")) %>% 
  select(group, feature, pathway, coef, pval, star1, star2, star3) %>% arrange(desc(coef))

# biomarker
pwy_bmi <- results_bmi %>% 
  mutate(group="BMI", 
         star1=case_when(fdr<0.05 ~ "+"),
         star2=case_when(fdr>=0.05 &  fdr<0.1 ~ "#"),
         star3=case_when(fdr>=0.1 &  fdr<0.25 ~ "*")) %>% 
  select(group, feature, pathway, coef, pval, star1, star2, star3)
pwy_ldlc <- results_ldlc %>% 
  mutate(group="LDL-C", 
         star1=case_when(fdr<0.05 ~ "+"),
         star2=case_when(fdr>=0.05 &  fdr<0.1 ~ "#"),
         star3=case_when(fdr>=0.1 &  fdr<0.25 ~ "*")) %>% 
  select(group, feature, pathway, coef, pval, star1, star2, star3)
pwy_hdlc <- results_hdlc %>%
  mutate(group="HDL-C", 
         star1=case_when(fdr<0.05 ~ "+"),
         star2=case_when(fdr>=0.05 &  fdr<0.1 ~ "#"),
         star3=case_when(fdr>=0.1 &  fdr<0.25 ~ "*")) %>% 
  select(group, feature, pathway, coef, pval, star1, star2, star3)
pwy_tg <- results_tg %>%
  mutate(group="TG", 
         star1=case_when(fdr<0.05 ~ "+"),
         star2=case_when(fdr>=0.05 &  fdr<0.1 ~ "#"),
         star3=case_when(fdr>=0.1 &  fdr<0.25 ~ "*")) %>% 
  select(group, feature, pathway, coef, pval, star1, star2, star3)
pwy_homair <- results_homair %>% 
  mutate(group="HOMA-IR", 
         star1=case_when(fdr<0.05 ~ "+"),
         star2=case_when(fdr>=0.05 &  fdr<0.1 ~ "#"),
         star3=case_when(fdr>=0.1 &  fdr<0.25 ~ "*")) %>% 
  select(group, feature, pathway, coef, pval, star1, star2, star3)
pwy_homab <- results_homab %>%
  mutate(group="HOMA-B", 
         star1=case_when(fdr<0.05 ~ "+"),
         star2=case_when(fdr>=0.05 &  fdr<0.1 ~ "#"),
         star3=case_when(fdr>=0.1 &  fdr<0.25 ~ "*")) %>% 
  select(group, feature, pathway, coef, pval, star1, star2, star3)
pwy_crp <- results_crp %>% 
  mutate(group="hs-CRP", 
         star1=case_when(fdr<0.05 ~ "+"),
         star2=case_when(fdr>=0.05 &  fdr<0.1 ~ "#"),
         star3=case_when(fdr>=0.1 &  fdr<0.25 ~ "*")) %>% 
  select(group, feature, pathway, coef, pval, star1, star2, star3)

heatdata0 <- pwy_T2D

heatdata1 <- rbind(pwy_bmi, pwy_ldlc)
heatdata1 <- rbind(heatdata1, pwy_hdlc)
heatdata1 <- rbind(heatdata1, pwy_tg)
heatdata1 <- rbind(heatdata1, pwy_homair)
heatdata1 <- rbind(heatdata1, pwy_homab)
heatdata1 <- rbind(heatdata1, pwy_crp)

heatdata1$group <- factor(heatdata1$group, levels=c("BMI", "HOMA-IR", "HOMA-B", "HDL-C",
                                                    "LDL-C", "TG", "hs-CRP"))
summary(heatdata0$coef)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.06027  0.11762  0.15482  0.13667  0.18538  0.37970
summary(heatdata1$coef)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.169700 -0.016894  0.006611  0.007352  0.033509  0.169047 

sum(heatdata1$coef>0.15) # 1
sum(heatdata1$coef<(-0.15)) # 1
heatdata1$coef[heatdata1$coef>=0.15] <- 0.15
heatdata1$coef[heatdata1$coef<(-0.15)] <- (-0.15)

sum(pwy_T2D$coef>0) # 14
sum(pwy_T2D$coef<0) # 3
pwy1 <- pwy_T2D[pwy_T2D$coef>0, ] %>% arrange(desc(coef)) %>% pull(pathway)
pwy2 <- pwy_T2D[pwy_T2D$coef<0, ] %>% arrange(coef) %>% pull(pathway)

# cluster pwys
library(ComplexHeatmap)
library(circlize)

# positive
h_all <- heatdata1[,c("group", "pathway", "coef")]
hall <- reshape(h_all, idvar = "pathway", timevar = "group", direction = "wide")
rownames(hall) <- hall$pathway
hall <- hall[,-1]
colnames(hall) <- gsub("coef.", "", colnames(hall))

col_fun = colorRamp2(c(-0.4, 0, 0.4), c("#08519c", "white", "#bd0026"))
heatall <- Heatmap(as.matrix(hall), 
                 col = col_fun,
                 cluster_columns = FALSE)
heatall <- draw(heatall)
pwy_all <- rownames(hall)[row_order(heatall)]
pwy_all_dend <- row_dend(heatall)
pwy_all_dend_data <- dendro_data(pwy_all_dend, type = "rectangle")

# "LACTOSECAT-PWY: lactose and galactose degradation I" -> "LACTOSECAT-PWY: lactose degradation I"
# "NAGLIPASYN-PWY: lipid IVA biosynthesis" -> "NAGLIPASYN-PWY: lipid IVA biosynthesis (E.coli)"
# "PWY-5971: palmitate biosynthesis II (bacteria and plants)" -> "PWY-5971: palmitate biosynthesis II (type II fatty acid synthase)"
# "PWY-7357: thiamin formation from pyrithiamine and oxythiamine (yeast)" -> "PWY-7357: thiamine diphosphate formation from pyrithiamine and oxythiamine (yeast)"
# "TEICHOICACID-PWY: teichoic acid (poly-glycerol) biosynthesis" -> "TEICHOICACID-PWY: poly(glycerol phosphate) wall teichoic acid biosynthesis"

heatdata1$pathway[heatdata1$pathway=="LACTOSECAT-PWY: lactose and galactose degradation I"] <- "LACTOSECAT-PWY: lactose degradation I"
heatdata1$pathway[heatdata1$pathway=="NAGLIPASYN-PWY: lipid IVA biosynthesis"] <- "NAGLIPASYN-PWY: lipid IVA biosynthesis (E.coli)"
heatdata1$pathway[heatdata1$pathway=="PWY-5971: palmitate biosynthesis II (bacteria and plants)"] <- "PWY-5971: palmitate biosynthesis II (type II fatty acid synthase)"
heatdata1$pathway[heatdata1$pathway=="PWY-7357: thiamin formation from pyrithiamine and oxythiamine (yeast)"] <- "PWY-7357: thiamine diphosphate formation from pyrithiamine and oxythiamine (yeast)"
heatdata1$pathway[heatdata1$pathway=="TEICHOICACID-PWY: teichoic acid (poly-glycerol) biosynthesis"] <- "TEICHOICACID-PWY: poly(glycerol phosphate) wall teichoic acid biosynthesis"

pwy_all[pwy_all=="LACTOSECAT-PWY: lactose and galactose degradation I"] <- "LACTOSECAT-PWY: lactose degradation I"
pwy_all[pwy_all=="NAGLIPASYN-PWY: lipid IVA biosynthesis"] <- "NAGLIPASYN-PWY: lipid IVA biosynthesis (E.coli)"
pwy_all[pwy_all=="PWY-5971: palmitate biosynthesis II (bacteria and plants)"] <- "PWY-5971: palmitate biosynthesis II (type II fatty acid synthase)"
pwy_all[pwy_all=="PWY-7357: thiamin formation from pyrithiamine and oxythiamine (yeast)"] <- "PWY-7357: thiamine diphosphate formation from pyrithiamine and oxythiamine (yeast)"
pwy_all[pwy_all=="TEICHOICACID-PWY: teichoic acid (poly-glycerol) biosynthesis"] <- "TEICHOICACID-PWY: poly(glycerol phosphate) wall teichoic acid biosynthesis"

heatdata1$lab <- factor(heatdata1$pathway, levels=rev(pwy_all))

library(ggplot2)
library(ggdendro)
library(scales)
heat1_all <- ggplot(heatdata1, aes(x=group, y=lab)) +
  geom_point(aes(fill=coef), size=11, shape=22) +
  scale_fill_gradient2(low="#08519c",
                       mid="white",
                       high="#bd0026",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.15,0,0.15),
                       limits=c(-0.15,0.15)) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  geom_text(aes(label=star1), color="black", size=5, nudge_y = 0.05) +
  geom_text(aes(label=star2), color="black", size=5, nudge_y = 0.0) +
  geom_text(aes(label=star3), color="black", size=7, nudge_y = -0.2) +
  labs(color=expression(beta *" coefficients")) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.8),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(angle = 60, vjust = 0, hjust = 0, color = "black", size = 12),
        axis.ticks.length=unit(0.07,"inch"),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, size = 1.5))

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R1/figS8/figure_s8a.pdf", 
    width = 9, height = 6.5, onefile = F) # Open a new pdf file
print(heat1_all)
dev.off() # Close the file

