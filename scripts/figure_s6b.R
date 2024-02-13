library(tidyverse)
library(MMUPHin)

setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R1/figS6")
load(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s6/figure_s6c_ordinal.RData")

######################################################
#        13 species confounded by metf
######################################################
metf_plot13 <- meta_metf_results[meta_metf_results$feature %in% spe13,] %>% arrange(coef)

# data for p0_1
species0 <- metf_plot13 %>% mutate(lab = gsub("s__", "", feature))
species0$lab <- gsub("_", " ", species0$lab)
species0$lab <- factor(species0$lab, levels=species0$lab)
species0$lci <- species0$coef-species0$stderr
species0$uci <- species0$coef+species0$stderr

# data for p0_2
species_0 <- metf_plot13$feature
karlsson0 <- meta_metf$maaslin_fits$`Karlsson_2013 (SWE)` %>% filter(feature %in% species_0) %>% mutate(study="Karlsson_2013 (SWE)")
direct0 <- meta_metf$maaslin_fits$`DIRECT-PLUS (ISR)` %>% filter(feature %in% species_0) %>% mutate(study="DIRECT-PLUS (ISR)")
mbs0 <- meta_metf$maaslin_fits$`NHSII (USA)` %>% filter(feature %in% species_0) %>% mutate(study="NHSII (USA)")
mlvs0 <- meta_metf$maaslin_fits$`HPFS (USA)` %>% filter(feature %in% species_0) %>% mutate(study="HPFS (USA)")
qin0 <- meta_metf$maaslin_fits$`Qin_2012 (CHN)` %>% filter(feature %in% species_0) %>% mutate(study="Qin_2012 (CHN)")
sol0 <- meta_metf$maaslin_fits$`SOL (USA)` %>% filter(feature %in% species_0) %>% mutate(study="SOL (USA)")
pedersen0 <- meta_metf$maaslin_fits$`Fromentin_2022 (DEU_DNK_FRA)` %>% filter(feature %in% species_0) %>% mutate(study="Fromentin_2022 (DEU/DNK/FRA)")

studies0 <- rbind(karlsson0, qin0)
studies0 <- rbind(studies0, mbs0)
studies0 <- rbind(studies0, mlvs0)
studies0 <- rbind(studies0, direct0)
studies0 <- rbind(studies0, sol0)
studies0 <- rbind(studies0, pedersen0)
studies0 <- studies0[!is.na(studies0$coef),]

studies0$study <- ifelse(studies0$study=="SOL (USA)",
                         "HCHS/SOL (USA)",
                         studies0$study)

studies0$study <- factor(studies0$study, levels=c("DIRECT-PLUS (ISR)","Fromentin_2022 (DEU/DNK/FRA)",
                                                  "HPFS (USA)","Karlsson_2013 (SWE)", 
                                                  "NHSII (USA)", "Qin_2012 (CHN)","HCHS/SOL (USA)"))
studies0 <- studies0 %>% mutate(star1=case_when(qval<0.05 ~ "+"),
                                star2=case_when(qval>=0.05 &  qval<0.1 ~ "#"),
                                star3=case_when(qval>=0.1 &  qval<0.25 ~ "*"))
studies0$feature <- factor(studies0$feature, levels=species0$feature)

# coef
summary(studies0$coef)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -5.7603 -0.9931 -0.5340 -0.5197  0.2704  2.2577
sum(studies0$coef<(-3)) # 5
sum(studies0$coef>3) # 0
studies0$coef[studies0$coef>=3] <- 3
studies0$coef[studies0$coef<=(-3)] <- (-3)

summary(species0$lci)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.3825 -1.0826 -0.6747 -0.4211  0.3391  1.1744
summary(species0$uci)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.84587 -0.70463 -0.29113  0.03399  0.65191  1.83468

species0_13 <- species0
studies0_13 <- studies0

library(ggplot2)
p13_1 <- ggplot(species0_13, aes(x=coef, y=lab, xmin=lci, xmax=uci)) +
  geom_errorbar(width=0.3, size=1) +
  geom_point(aes(color=fdr), size=3) +
  scale_color_gradient(low="#132B43", high="#56B1F7", 
                       limits=c(0,0.25), breaks=c(0, 0.25)) +
  geom_vline(xintercept=0, linetype="dashed", size=1.2) +
  scale_x_continuous(breaks=c(-1,0,1,2), limits=c(-1.5, 2.2)) +
  labs(x=expression(beta~coefficient~(SE)),
       color="FDR") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "italic", color = "black", size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top",
        legend.text = element_text(colour = "black", size = 16),
        legend.title = element_text(colour = "black", size = 16)) +
  guides(color = guide_colourbar(title.position="top", 
                                 title.hjust = 0.5,
                                 barheight = 1.5,
                                 barwidth = 10))

library(scales)
p13_2 <- ggplot(studies0_13, aes(x=study, y=feature)) +
  geom_point(aes(fill=coef), size=10, shape=21, stroke=0) +
  scale_fill_gradientn(
    colors=c("#08519c","white","#bd0026"),
    values=rescale(c(-3,0,3)),
    breaks=c(-3,0,3),
    limits=c(-3,3)) +
  geom_text(aes(label=star1), color="black", size=5, nudge_y = 0.05) +
  geom_text(aes(label=star2), color="black", size=5, nudge_y = 0.0) +
  geom_text(aes(label=star3), color="black", size=7, nudge_y = -0.2) +
  labs(fill="Coefficients for\n species-metf") +
  xlab("") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "top",
        legend.text = element_text(colour = "black", size = 16),
        legend.title = element_text(colour = "black", size = 16)) +
  guides(fill = guide_colourbar(title.position="top", 
                                title.hjust = 0.5,
                                barheight = 1.5,
                                barwidth = 10))

# panel 3
data1 <- meta_1[meta_1$feature %in% spe13, ]
data3 <- meta_3[meta_3$feature %in% spe13, ]
data1$model <- "Not adjusting for metf"
data3$model <- "Adjusting for metf"
heatmap <- rbind(data1[, c("feature", "model", "coef", "pval", "qval.fdr")],
                 data3[, c("feature", "model", "coef", "pval", "qval.fdr")])
heatmap <- heatmap %>% mutate(star1=case_when(qval.fdr<0.05 ~ "+"),
                              star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "#"),
                              star3=case_when(qval.fdr>=0.1 &  qval.fdr<0.25 ~ "*"))
heatmap$feature <- factor(heatmap$feature, levels=species0_13$feature)
summary(heatmap$coef)
#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# -0.33704 -0.17022 -0.04567 -0.04543  0.03378  0.37283 

sum(heatmap$coef>=0.5) # 0
sum(heatmap$coef<(-0.4)) # 0
# heatmap$coef[heatmap$coef>=0.5] <- 0.5
# heatmap$coef[heatmap$coef<(-0.4)] <- (-0.4)
heatmap$model <- factor(heatmap$model, levels=c("Not adjusting for metf",
                                                "Adjusting for metf"))

heatmap_13 <- heatmap
p13_3 <- ggplot(heatmap_13, aes(x=model, y=feature)) +
  geom_point(aes(color=coef), size=10, shape=15) +
  scale_color_gradient2(low="#08519c",
                        mid="white",
                        high="#bd0026",
                        midpoint = 0,
                        guide = "colourbar",
                        breaks=c(-0.4,0,0.5),
                        limits=c(-0.4,0.5)) +
  geom_text(aes(label=star1), color="black", size=5, nudge_y = 0.05) +
  geom_text(aes(label=star2), color="black", size=5, nudge_y = 0.0) +
  geom_text(aes(label=star3), color="black", size=7, nudge_y = -0.2) +
  labs(color="Coefficients for\n species-T2D") +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        legend.position = "right",
        legend.direction="horizontal",
        legend.text = element_text(colour = "black", size = 16),
        legend.title = element_text(colour = "black", size = 16)) +
  guides(color = guide_colourbar(title.position="top", 
                                 title.hjust = 0.5,
                                 barheight = 1.5,
                                 barwidth = 10))
#  theme(plot.margin = margin(0.1, 3, 0.1, 0.1, "cm"))

######################################################
#    7 species associated with both T2D and metf
######################################################
metf_plot7 <- meta_metf_results[meta_metf_results$feature %in% spe7,] %>% arrange(coef)

# data for p0_1
species0 <- metf_plot7 %>% mutate(lab = gsub("s__", "", feature))
species0$lab <- gsub("_", " ", species0$lab)
species0$lab <- factor(species0$lab, levels=species0$lab)
species0$lci <- species0$coef-species0$stderr
species0$uci <- species0$coef+species0$stderr

# data for p0_2
species_0 <- metf_plot7$feature
karlsson0 <- meta_metf$maaslin_fits$`Karlsson_2013 (SWE)` %>% filter(feature %in% species_0) %>% mutate(study="Karlsson_2013 (SWE)")
direct0 <- meta_metf$maaslin_fits$`DIRECT-PLUS (ISR)` %>% filter(feature %in% species_0) %>% mutate(study="DIRECT-PLUS (ISR)")
mbs0 <- meta_metf$maaslin_fits$`NHSII (USA)` %>% filter(feature %in% species_0) %>% mutate(study="NHSII (USA)")
mlvs0 <- meta_metf$maaslin_fits$`HPFS (USA)` %>% filter(feature %in% species_0) %>% mutate(study="HPFS (USA)")
qin0 <- meta_metf$maaslin_fits$`Qin_2012 (CHN)` %>% filter(feature %in% species_0) %>% mutate(study="Qin_2012 (CHN)")
sol0 <- meta_metf$maaslin_fits$`SOL (USA)` %>% filter(feature %in% species_0) %>% mutate(study="SOL (USA)")
pedersen0 <- meta_metf$maaslin_fits$`Fromentin_2022 (DEU_DNK_FRA)` %>% filter(feature %in% species_0) %>% mutate(study="Fromentin_2022 (DEU/DNK/FRA)")

studies0 <- rbind(karlsson0, qin0)
studies0 <- rbind(studies0, mbs0)
studies0 <- rbind(studies0, mlvs0)
studies0 <- rbind(studies0, direct0)
studies0 <- rbind(studies0, sol0)
studies0 <- rbind(studies0, pedersen0)
studies0 <- studies0[!is.na(studies0$coef),]

studies0$study[studies0$study=="SOL (USA)"] <- "HCHS/SOL (USA)"
studies0$study <- factor(studies0$study, levels=c("DIRECT-PLUS (ISR)","Fromentin_2022 (DEU/DNK/FRA)",
                                                  "HPFS (USA)","Karlsson_2013 (SWE)", 
                                                  "NHSII (USA)", "Qin_2012 (CHN)","HCHS/SOL (USA)"))
studies0 <- studies0 %>% mutate(star1=case_when(qval<0.05 ~ "+"),
                                star2=case_when(qval>=0.05 &  qval<0.1 ~ "#"),
                                star3=case_when(qval>=0.1 &  qval<0.25 ~ "*"))
studies0$feature <- factor(studies0$feature, levels=species0$feature)

# coef
summary(studies0$coef)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -2.91193 -0.79763 -0.35107 -0.08666  0.58840  4.59159 
sum(studies0$coef<(-3)) # 0
sum(studies0$coef>3) # 2
studies0$coef[studies0$coef>=3] <- 3
studies0$coef[studies0$coef<=(-3)] <- (-3)

summary(species0$lci)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.79337 -0.66678 -0.61111  0.09159  0.84378  1.61040
summary(species0$uci)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.3270 -0.2697 -0.2359  0.5597  1.5477  2.0832

species0_7 <- species0
studies0_7 <- studies0
p7_1 <- ggplot(species0_7, aes(x=coef, y=lab, xmin=lci, xmax=uci)) +
  geom_errorbar(width=0.3, size=1) +
  geom_point(aes(color=fdr), size=3) +
  scale_color_gradient(low="#132B43", high="#56B1F7", 
                       limits=c(0,0.25), breaks=c(0, 0.25)) +
  geom_vline(xintercept=0, linetype="dashed", size=1.2) +
  scale_x_continuous(breaks=c(-1,0,1,2), limits=c(-1.5, 2.2)) +
  labs(x=expression(beta~coefficient~(SE)),
       color="FDR") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_text(face = "italic", color = "black", size = 16),
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_blank(),
        legend.position = "none")

p7_2 <- ggplot(studies0, aes(x=study, y=feature)) +
  geom_point(aes(fill=coef), size=10, shape=21, stroke=0) +
  scale_fill_gradientn(
    colors=c("#08519c","white","#bd0026"),
    values=rescale(c(-3,0,3)),
    breaks=c(-3,0,3),
    limits=c(-3,3)) +
  geom_text(aes(label=star1), color="black", size=5, nudge_y = 0.05) +
  geom_text(aes(label=star2), color="black", size=5, nudge_y = 0.0) +
  geom_text(aes(label=star3), color="black", size=7, nudge_y = -0.2) +
  labs(fill=expression(beta * " coefficient")) +
  xlab("") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1, size=16, colour = "black"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")

# panel 3
data1 <- meta_1[meta_1$feature %in% spe7, ]
data3 <- meta_3[meta_3$feature %in% spe7, ]
data1$model <- "Not adjusting for metf"
data3$model <- "Adjusting for metf"
heatmap <- rbind(data1[, c("feature", "model", "coef", "pval", "qval.fdr")],
                 data3[, c("feature", "model", "coef", "pval", "qval.fdr")])
heatmap <- heatmap %>% mutate(star1=case_when(qval.fdr<0.05 ~ "+"),
                              star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "#"),
                              star3=case_when(qval.fdr>=0.1 &  qval.fdr<0.25 ~ "*"))
heatmap$feature <- factor(heatmap$feature, levels=species0_7$feature)
summary(heatmap$coef)
#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# -0.47354 -0.30545 -0.22002 -0.06206  0.18063  0.71066

sum(heatmap$coef>=0.5) # 2
sum(heatmap$coef<(-0.4)) # 1
heatmap$coef[heatmap$coef>=0.5] <- 0.5
heatmap$coef[heatmap$coef<(-0.4)] <- (-0.4)
heatmap$model <- factor(heatmap$model, levels=c("Not adjusting for metf",
                                                "Adjusting for metf"))

heatmap_7 <- heatmap
p7_3 <- ggplot(heatmap_7, aes(x=model, y=feature)) +
  geom_point(aes(color=coef), size=10, shape=15) +
  scale_color_gradient2(low="#08519c",
                        mid="white",
                        high="#bd0026",
                        midpoint = 0,
                        guide = "colourbar",
                        breaks=c(-0.4,0,0.5),
                        limits=c(-0.4,0.5)) +
  geom_text(aes(label=star1), color="black", size=5, nudge_y = 0.05) +
  geom_text(aes(label=star2), color="black", size=5, nudge_y = 0.0) +
  geom_text(aes(label=star3), color="black", size=7, nudge_y = -0.2) +
  labs(color=expression(beta * " coefficient")) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.6),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1, size=16, colour = "black"),
        axis.ticks.length=unit(0.07,"inch"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5)) +
  theme(legend.position = "none")

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R1/figS6/figure_s6b.pdf", 
    width = 12, height = 11.5, onefile = F) # Open a new pdf file
egg::ggarrange(p13_1, p13_2, p13_3, 
               p7_1, p7_2, p7_3, 
               nrow=2, ncol=3, 
               widths = c(2,1.6,0.5),
               heights = c(13,7))
dev.off() # Close the file
