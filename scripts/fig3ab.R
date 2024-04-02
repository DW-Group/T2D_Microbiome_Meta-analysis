setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig3")

#########################################################
#                     pathways
#########################################################
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/analysis/abund0020/pathway_trend.RData")
# select pwys to plot
meta_3_plot <- meta_3trend %>% 
  filter(feature %in% c(
    "GLYCOLYSIS",
    "HEXITOLDEGSUPER-PWY",
    "METHGLYUT-PWY",
    "PENTOSE-P-PWY",
    "SER-GLYSYN-PWY",
    "FASYN-ELONG-PWY",
    "PWY-5347",
    "PWY-5971",
    "PWY-6901",
    "PWY-7357",
    "SULFATE-CYS-PWY",
    "TEICHOICACID-PWY",
    "NAGLIPASYN-PWY"
  )) %>% mutate(group=case_when(feature=="PWY-7357" ~ "Vitamin biosynthesis",
                                feature=="PWY-6901" ~ "Glucose homeostasis",
                                feature=="PWY-5971" ~ "Saturated fatty acid biosynthesis",
                                feature=="FASYN-ELONG-PWY" ~ "Saturated fatty acid biosynthesis",
                                feature=="SER-GLYSYN-PWY" ~ "Sulfur-containing amino acid biosynthesis",
                                feature=="SULFATE-CYS-PWY" ~ "Sulfur-containing amino acid biosynthesis",
                                feature=="TEICHOICACID-PWY" ~ "Biosynthesis of bacterial structural components",
                                feature=="GLYCOLYSIS" ~ "Glucose homeostasis",
                                feature=="HEXITOLDEGSUPER-PWY" ~ "Glucose homeostasis",
                                feature=="METHGLYUT-PWY" ~ "Glucose homeostasis",
                                feature=="NAGLIPASYN-PWY" ~ "Biosynthesis of bacterial structural components",
                                feature=="PENTOSE-P-PWY" ~ "Glucose homeostasis",
                                feature=="PWY-5347" ~ "Sulfur-containing amino acid biosynthesis"))
meta_3_plot$group <- factor(meta_3_plot$group)
meta_3_plot <- meta_3_plot %>% arrange(desc(group), coef)
meta_3_plot$pathway[meta_3_plot$pathway=="LACTOSECAT-PWY: lactose and galactose degradation I"] <- "LACTOSECAT-PWY: lactose degradation I"
meta_3_plot$pathway[meta_3_plot$pathway=="NAGLIPASYN-PWY: lipid IVA biosynthesis"] <- "NAGLIPASYN-PWY: lipid IVA biosynthesis (E.coli)"
meta_3_plot$pathway[meta_3_plot$pathway=="PWY-5971: palmitate biosynthesis II (bacteria and plants)"] <- "PWY-5971: palmitate biosynthesis II (type II fatty acid synthase)"
meta_3_plot$pathway[meta_3_plot$pathway=="PWY-7357: thiamin formation from pyrithiamine and oxythiamine (yeast)"] <- "PWY-7357: thiamine diphosphate formation from pyrithiamine and oxythiamine (yeast)"
meta_3_plot$pathway[meta_3_plot$pathway=="TEICHOICACID-PWY: teichoic acid (poly-glycerol) biosynthesis"] <- "TEICHOICACID-PWY: poly(glycerol phosphate) wall teichoic acid biosynthesis"

# data for p0_1
pathway0 <- meta_3_plot %>% mutate(lab = pathway)
pathway0$lab <- factor(pathway0$lab, levels=pathway0$lab)
pathway0$lci <- pathway0$coef-pathway0$stderr
pathway0$uci <- pathway0$coef+pathway0$stderr

# data for p0_2
pathway_0 <- meta_3_plot$feature
karlsson0 <- meta3_trend$maaslin_fits$`Karlsson_2013 (SWE)` %>% filter(feature %in% pathway_0) %>% mutate(study="Karlsson_2013 (SWE)")
direct0 <- meta3_trend$maaslin_fits$`DIRECT-PLUS (ISR)` %>% filter(feature %in% pathway_0) %>% mutate(study="DIRECT-PLUS (ISR)")
mbs0 <- meta3_trend$maaslin_fits$`NHSII (USA)` %>% filter(feature %in% pathway_0) %>% mutate(study="NHSII (USA)")
mlvs0 <- meta3_trend$maaslin_fits$`HPFS (USA)` %>% filter(feature %in% pathway_0) %>% mutate(study="HPFS (USA)")
qin0 <- meta3_trend$maaslin_fits$`Qin_2012 (CHN)` %>% filter(feature %in% pathway_0) %>% mutate(study="Qin_2012 (CHN)")
zhong0 <- meta3_trend$maaslin_fits$`Zhong_2019 (CHN)` %>% filter(feature %in% pathway_0) %>% mutate(study="Zhong_2019 (CHN)")
sol0 <- meta3_trend$maaslin_fits$`SOL (USA)` %>% filter(feature %in% pathway_0) %>% mutate(study="SOL (USA)")
pedersen0 <- meta3_trend$maaslin_fits$`Fromentin_2022 (DEU_DNK_FRA)` %>% filter(feature %in% pathway_0) %>% mutate(study="Fromentin_2022 (DEU/DNK/FRA)")
wua0 <- meta3_trend$maaslin_fits$`Wu_2020a (SWE)` %>% filter(feature %in% pathway_0) %>% mutate(study="Wu_2020a (SWE)")
wub0 <- meta3_trend$maaslin_fits$`Wu_2020b (SWE)` %>% filter(feature %in% pathway_0) %>% mutate(study="Wu_2020b (SWE)")

karlsson0 <- merge(pathway0[, c("feature", "pathway")], 
                   karlsson0, by="feature", all.y=T)
direct0 <- merge(pathway0[, c("feature", "pathway")], 
                 direct0, by="feature", all.y=T)
mbs0 <- merge(pathway0[, c("feature", "pathway")], 
              mbs0, by="feature", all.y=T)
mlvs0 <- merge(pathway0[, c("feature", "pathway")], 
               mlvs0, by="feature", all.y=T)
qin0 <- merge(pathway0[, c("feature", "pathway")], 
              qin0, by="feature", all.y=T)
zhong0 <- merge(pathway0[, c("feature", "pathway")], 
                zhong0, by="feature", all.y=T)
sol0 <- merge(pathway0[, c("feature", "pathway")], 
              sol0, by="feature", all.y=T)
pedersen0 <- merge(pathway0[, c("feature", "pathway")], 
                   pedersen0, by="feature", all.y=T)
wua0 <- merge(pathway0[, c("feature", "pathway")], 
             wua0, by="feature", all.y=T)
wub0 <- merge(pathway0[, c("feature", "pathway")], 
              wub0, by="feature", all.y=T)

studies0 <- rbind(karlsson0, qin0)
studies0 <- rbind(studies0, mbs0)
studies0 <- rbind(studies0, mlvs0)
studies0 <- rbind(studies0, zhong0)
studies0 <- rbind(studies0, direct0)
studies0 <- rbind(studies0, sol0)
studies0 <- rbind(studies0, pedersen0)
studies0 <- rbind(studies0, wua0)
studies0 <- rbind(studies0, wub0)
studies0 <- studies0[!is.na(studies0$coef),]

studies0$study <- ifelse(studies0$study=="SOL (USA)",
                         "HCHS/SOL (USA)",
                         studies0$study)
studies0$study <- factor(studies0$study, levels=c("DIRECT-PLUS (ISR)", "Fromentin_2022 (DEU/DNK/FRA)", "HPFS (USA)", 
                                                  "Karlsson_2013 (SWE)", "NHSII (USA)", "Qin_2012 (CHN)",
                                                  "HCHS/SOL (USA)", "Wu_2020a (SWE)", "Wu_2020b (SWE)", "Zhong_2019 (CHN)"))
studies0 <- studies0 %>% mutate(star1=case_when(qval<0.05 ~ "#"),
                                star2=case_when(qval>=0.05 &  qval<0.1 ~ "*"))
studies0$pathway <- factor(studies0$pathway, levels=pathway0$pathway)

# coef
summary(studies0$coef)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.47027 -0.01874  0.06991  0.07047  0.15258  0.69705
sum(studies0$coef<(-0.3)) # 3
sum(studies0$coef>0.6) # 2
studies0$coef[studies0$coef>=0.6] <- 0.6
studies0$coef[studies0$coef<=(-0.3)] <- (-0.3)

summary(pathway0$lci)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.03460  0.03843  0.04297  0.05191  0.06192  0.13516
summary(pathway0$uci)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.01143  0.08363  0.09564  0.11134  0.14549  0.21204

pwy_1_data <- pathway0
pwy_2_data <- studies0

library(ggplot2)
pwy_1_data$coef2 <- sign(pwy_1_data$coef)*sqrt(abs(pwy_1_data$coef))
pwy_1_data$lci2 <- sign(pwy_1_data$lci)*sqrt(abs(pwy_1_data$lci))
pwy_1_data$uci2 <- sign(pwy_1_data$uci)*sqrt(abs(pwy_1_data$uci))
brk <- c(-0.05,0,0.1,0.2)
brk <- sign(brk)*sqrt(abs(brk))
pwy_1 <- ggplot(pwy_1_data, aes(x=coef2, y=lab, xmin=lci2, xmax=uci2)) +
  geom_errorbar(width=0.3, size=1) +
  geom_point(aes(color=qval.fdr), size=3) +
  scale_color_gradient(low="#56B1F7", high="#132B43", 
                       limits=c(0,0.1), breaks=c(0, 0.1)) +
  geom_vline(xintercept=0, linetype="dashed", size=1.2) +
  scale_x_continuous(breaks=brk,
                     labels=c(-0.05,0,0.1,0.2), 
                     limits=c(brk[1],0.50)) +
  labs(x="",
       color="FDR") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1.2),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 16),
        axis.title = element_blank(),
        axis.title.y=element_blank(),
        legend.position = c(0.2,0.79),
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14))

library(scales)
pwy_2 <- ggplot(pwy_2_data, aes(x=study, y=pathway)) +
  geom_point(aes(fill=coef), size=10, shape=21, stroke=0.5) +
  scale_fill_gradientn(
    colors=c("#08519c","white","#bd0026"),
    values=rescale(c(-0.3,0,0.6)),
    breaks=c(-0.3,0,0.6),
    limits=c(-0.3,0.6)) +
  geom_text(aes(label=star1), color="black", size=7, nudge_y = 0.0) +
  geom_text(aes(label=star2), color="black", size=9, nudge_y = -0.2) +
  labs(fill=expression(beta * " coefficient (pathways)")) +
  xlab("") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1.2),
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
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig4/pwy_spe_v3.RData")

length(unique(pwy_all_spe_v3_xx$species)) # 15
length(unique(pwy_all_spe_v3_xx$path)) # 17
summary(pwy_all_spe_v3_xx$avg)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -4.7681 -0.0699  1.2139  1.2610  3.0188  9.4003

pwy_all_spe_v3_xx$pathway <- factor(pwy_all_spe_v3_xx$pathway,
                                    levels = pathway0$pathway)
pwy_all_spe_v3_xx$species <- gsub("s__", "", pwy_all_spe_v3_xx$species)
pwy_all_spe_v3_xx$species <- gsub("_", " ", pwy_all_spe_v3_xx$species)

gdata::keep(pwy_1_data, pwy_2_data, brk,
            pwy_1, pwy_2, pwy_all_spe_v3_xx, sure=T)

#########################################################
#                     enzymes
#########################################################
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/analysis/abund0020/enzyme_trend.RData")
# select ecs to plot
meta_3_plot <- meta_3trend %>% filter(qval.fdr<0.1) %>% 
  filter(feature %in% c(
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
    "EC2.4.1.345: NO_NAME",
    "EC3.4.24.55: Pitrilysin",
    "EC1.2.1.12: Glyceraldehyde-3-phosphate dehydrogenase (phosphorylating)",
    "EC1.1.1.136: UDP-N-acetylglucosamine 6-dehydrogenase",
    "EC1.1.1.44: Phosphogluconate dehydrogenase (NADP(+)-dependent, decarboxylating)",
    "EC1.1.1.49: Glucose-6-phosphate dehydrogenase (NADP(+))",
    "EC3.1.6.6: Choline-sulfatase",
    "EC2.7.8.31: Undecaprenyl-phosphate glucose phosphotransferase"
  )) %>%  mutate(group=case_when(feature %in% c("EC2.4.2.43: Lipid IV(A) 4-amino-4-deoxy-L-arabinosyltransferase",
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
meta_3_plot$group <- factor(meta_3_plot$group)
meta_3_plot <- meta_3_plot %>% arrange(desc(group), coef)

# data for p0_1
enzyme0 <- meta_3_plot %>%  mutate(lab = feature)
enzyme0$lab[enzyme0$lab=="EC2.4.1.345: NO_NAME"] <- "EC2.4.1.345: Phosphatidyl-myo-inositol a-mannosyltransferase"
enzyme0$lab <- factor(enzyme0$lab, levels=enzyme0$lab)
enzyme0$lci <- enzyme0$coef-enzyme0$stderr
enzyme0$uci <- enzyme0$coef+enzyme0$stderr

# data for p0_2
enzyme_0 <- meta_3_plot$feature
karlsson0 <- meta3_trend$maaslin_fits$`Karlsson_2013 (SWE)` %>% filter(feature %in% enzyme_0) %>% mutate(study="Karlsson_2013 (SWE)")
direct0 <- meta3_trend$maaslin_fits$`DIRECT-PLUS (ISR)` %>% filter(feature %in% enzyme_0) %>% mutate(study="DIRECT-PLUS (ISR)")
mbs0 <- meta3_trend$maaslin_fits$`NHSII (USA)` %>% filter(feature %in% enzyme_0) %>% mutate(study="NHSII (USA)")
mlvs0 <- meta3_trend$maaslin_fits$`HPFS (USA)` %>% filter(feature %in% enzyme_0) %>% mutate(study="HPFS (USA)")
qin0 <- meta3_trend$maaslin_fits$`Qin_2012 (CHN)` %>% filter(feature %in% enzyme_0) %>% mutate(study="Qin_2012 (CHN)")
zhong0 <- meta3_trend$maaslin_fits$`Zhong_2019 (CHN)` %>% filter(feature %in% enzyme_0) %>% mutate(study="Zhong_2019 (CHN)")
sol0 <- meta3_trend$maaslin_fits$`SOL (USA)` %>% filter(feature %in% enzyme_0) %>% mutate(study="SOL (USA)")
pedersen0 <- meta3_trend$maaslin_fits$`Fromentin_2022 (DEU_DNK_FRA)` %>% filter(feature %in% enzyme_0) %>% mutate(study="Fromentin_2022 (DEU/DNK/FRA)")
wua0 <- meta3_trend$maaslin_fits$`Wu_2020a (SWE)` %>% filter(feature %in% enzyme_0) %>% mutate(study="Wu_2020a (SWE)")
wub0 <- meta3_trend$maaslin_fits$`Wu_2020b (SWE)` %>% filter(feature %in% enzyme_0) %>% mutate(study="Wu_2020b (SWE)")

studies0 <- rbind(karlsson0, qin0)
studies0 <- rbind(studies0, mbs0)
studies0 <- rbind(studies0, mlvs0)
studies0 <- rbind(studies0, zhong0)
studies0 <- rbind(studies0, direct0)
studies0 <- rbind(studies0, sol0)
studies0 <- rbind(studies0, pedersen0)
studies0 <- rbind(studies0, wua0)
studies0 <- rbind(studies0, wub0)
studies0 <- studies0[!is.na(studies0$coef),]

studies0$study <- ifelse(studies0$study=="SOL (USA)",
                         "HCHS/SOL (USA)",
                         studies0$study)
studies0$study <- factor(studies0$study, levels=c("DIRECT-PLUS (ISR)", "Fromentin_2022 (DEU/DNK/FRA)", "HPFS (USA)", 
                                                  "Karlsson_2013 (SWE)", "NHSII (USA)", "Qin_2012 (CHN)",
                                                  "HCHS/SOL (USA)", "Wu_2020a (SWE)", "Wu_2020b (SWE)", "Zhong_2019 (CHN)"))
studies0 <- studies0 %>% mutate(star1=case_when(qval<0.05 ~ "#"),
                                star2=case_when(qval>=0.05 &  qval<0.1 ~ "*"))
studies0$feature <- factor(studies0$feature, levels=enzyme0$feature)

# coef
summary(studies0$coef)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.429730 -0.009586  0.061062  0.070373  0.135086  0.399535 
sum(studies0$coef<(-0.3)) # 1
sum(studies0$coef>0.4) # 2
studies0$coef[studies0$coef>=0.4] <- 0.4
studies0$coef[studies0$coef<=(-0.3)] <- (-0.3)

summary(enzyme0$lci)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.03213  0.04341  0.06265  0.05633  0.07879  0.14796 
summary(enzyme0$uci)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.01589  0.09031  0.12601  0.11289  0.15055  0.23112

ecs_1_data <- enzyme0
ecs_2_data <- studies0

ecs_1_data$coef2 <- sign(ecs_1_data$coef)*sqrt(abs(ecs_1_data$coef))
ecs_1_data$lci2 <- sign(ecs_1_data$lci)*sqrt(abs(ecs_1_data$lci))
ecs_1_data$uci2 <- sign(ecs_1_data$uci)*sqrt(abs(ecs_1_data$uci))
brk2 <- c(-0.05,0,0.2,0.4)
brk2 <- sign(brk2)*sqrt(abs(brk2))

ecs_1_data$lab2 <- gsub("EC", "EC ", ecs_1_data$lab)
ecs_1 <- ggplot(ecs_1_data, aes(x=coef2, y=lab, xmin=lci2, xmax=uci2)) +
  geom_errorbar(width=0.3, size=1) +
  geom_point(aes(color=qval.fdr), size=3) +
  scale_color_gradient(low="#56B1F7", high="#132B43", 
                       limits=c(0,0.1), breaks=c(0, 0.1)) +
  geom_vline(xintercept=0, linetype="dashed", size=1.2) +
  scale_y_discrete(labels=ecs_1_data$lab2) +
  scale_x_continuous(breaks=brk2,
                     labels=c(-0.05,0,0.2,0.4), 
                     limits=c(brk2[1],0.68)) +
  labs(x=expression(beta~coefficient~(SE)),
       color="FDR") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1.2),
        axis.text.y = element_blank(), 
        axis.text.x = element_text(color = "black", size = 16),
        axis.title.x = element_text(colour = "black", hjust=0.5, vjust=(-1), size = 16),
        axis.title.y=element_blank(),
        legend.position = c(0.72,0.18),
        legend.justification = "left",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14)) 

ecs_2 <- ggplot(ecs_2_data, aes(x=study, y=feature)) +
  geom_point(aes(fill=coef), size=10, shape=21, stroke=0.5) +
  scale_fill_gradientn(
    colors=c("#08519c","white","#bd0026"),
    values=rescale(c(-0.3,0,0.4)),
    breaks=c(-0.3,0,0.4),
    limits=c(-0.3,0.4)) +
  geom_text(aes(label=star1), color="black", size=7, nudge_y = 0.0) +
  geom_text(aes(label=star2), color="black", size=9, nudge_y = -0.2) +
  labs(fill=expression(beta * " coefficient (enzymes)")) +
  xlab("") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1.2),
        axis.text.x = element_text(angle = 55, vjust = 1, hjust = 1, size=16, colour = "black"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "bottom",
        legend.text = element_text(colour = "black", size = 16),
        legend.title = element_text(colour = "black", size = 16)) +
  guides(fill = guide_colourbar(title.position="top", 
                                title.hjust = 0.5,
                                barheight = 1.5,
                                barwidth = 10))

# panel 3
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig4/ecs_spe_v3.RData")

length(unique(ecs_all_spe_v3_xx$species)) # 22
length(unique(ecs_all_spe_v3_xx$enzy)) # 19
summary(ecs_all_spe_v3_xx$avg)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -3.9195 -1.2099  0.8931  0.3336  1.3712  5.3409

ecs_all_spe_v3_xx$enzy <- factor(ecs_all_spe_v3_xx$enzy,
                                 levels = enzyme0$feature)
ecs_all_spe_v3_xx$species <- gsub("s__", "", ecs_all_spe_v3_xx$species)
ecs_all_spe_v3_xx$species <- gsub("_", " ", ecs_all_spe_v3_xx$species)

gdata::keep(pwy_1_data, pwy_2_data, brk,
            pwy_1, pwy_2, pwy_all_spe_v3_xx, ecs_1_data, ecs_2_data, brk2,
            ecs_1, ecs_2, ecs_all_spe_v3_xx, sure=T)

species <- sort(union(unique(pwy_all_spe_v3_xx$species),
                      unique(ecs_all_spe_v3_xx$species)))

library(paletteer) # https://r-charts.com/color-palettes/#discrete
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig4/figure4_col27.RData")
col <- col27$col
names(col) <- col27$species

pwy_3 <- ggplot(pwy_all_spe_v3_xx, aes(x=pathway, y=avg)) + 
  geom_hline(yintercept = 0, linetype="dashed", size=1.2) +
  geom_col(aes(fill=species), position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=col) +
  scale_y_continuous(breaks=c(-10,-5,0,5,10), limits=c(-11,11)) +
  guides(fill=guide_legend(ncol =1)) +
  coord_flip() +
  labs(y="") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1.2),
        axis.text.x = element_text(size=16, colour = "black"),
        axis.title.x = element_text(hjust=0.5, size=16, colour = "black"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_blank())

ecs_3 <- ggplot(ecs_all_spe_v3_xx, aes(x=enzy, y=avg)) + 
  geom_hline(yintercept = 0, linetype="dashed", size=1.2) +
  geom_col(aes(fill=species), position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=col) +
  scale_y_continuous(breaks=c(-4,0,4,8), limits=c(-6,9)) +
  guides(fill=guide_legend(ncol =1)) +
  coord_flip() +
  labs(y=expression("\n\n Mean difference\n in contribution %")) +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1.2),
        axis.text.x = element_text(size=16, colour = "black"),
        axis.title.x = element_text(hjust=0.5, vjust=(-21), size=16, colour = "black"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_blank())

colgroup <- c("#CC9900FF","#A9A9A9FF","#BA6338FF","#CE3D32FF",
              "#5DB1DDFF","#7A65A5FF","#339900FF","#99CC00FF")
names(colgroup) <- c("Biosynthesis of bacterial structural components", 
                     "Catabolism of dairy components", "Choline production", "Glucose homeostasis", 
                     "Saturated fatty acid biosynthesis", "Sulfur-containing amino acid biosynthesis",
                     "Tryptophan biosynthesis", "Vitamin biosynthesis")

pwy_0 <- ggplot(pwy_1_data, aes(x=1, y=lab)) +
  geom_point(aes(color=group), size=11.3, shape=15) +
  scale_color_manual(values=colgroup) +
  xlab("") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        # axis.line.y = element_line(size=1.2),
        axis.ticks.y = element_line(colour = "black", size = 1.2),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 17),
        axis.text.x = element_blank(),
        axis.title =element_blank(),
        legend.position = "none")

ecs_0 <- ggplot(ecs_1_data, aes(x=1, y=lab)) +
  geom_point(aes(color=group), size=11.3, shape=15) +
  labs(color="Functional pathways") +
  scale_color_manual(values=colgroup) +
  xlab("") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        # axis.line.y = element_line(size=1.2),
        axis.ticks.y = element_line(colour = "black", size = 1.2),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 17),
        axis.text.x = element_blank(),
        axis.title =element_blank(),
        legend.position = "none")

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig3/fig3ab.pdf", 
    width = 19.5, height = 16.4, onefile = F) # Open a new pdf file
egg::ggarrange(pwy_0, pwy_1, pwy_2, pwy_3, 
               ecs_0, ecs_1, ecs_2, ecs_3, 
               nrow=2, ncol=4,
               widths = c(0.06,1,1,1), heights = c(1.55,2))
dev.off() # Close the file

ecs_1_data$group <- as.character(ecs_1_data$group)
ecs_1_data[17:19, "group"] <- c("Saturated fatty acid biosynthesis", "Sulfur-containing amino acid biosynthesis",
                                "Vitamin biosynthesis")
ecs_1_data$group <- factor(ecs_1_data$group, 
                           levels=c("Biosynthesis of bacterial structural components", 
                                    "Choline production", "Glucose homeostasis", 
                                    "Saturated fatty acid biosynthesis", "Sulfur-containing amino acid biosynthesis",
                                    "Tryptophan biosynthesis", "Vitamin biosynthesis"))
# legend_grp
legend_grp <- ggplot(ecs_1_data, aes(x=1, y=lab)) +
  geom_point(aes(color=group), size=11.3, shape=15) +
  labs(color="Functional categories") +
  scale_color_manual(values=colgroup) +
  guides(color=guide_legend(ncol = 2, size = 2, byrow = T)) +
  xlab("") +
  theme_bw() +
  theme(axis.text = element_blank(),
        legend.position = "right",
        legend.text = element_text(colour = "black", size = 20),
        legend.title = element_text(colour = "black", size = 20),
        legend.key.height= unit(0.3, 'cm'),
        legend.key.width= unit(0.3, 'cm'),
        legend.spacing.x = unit(0.01, 'cm'),
        legend.spacing.y = unit(0.01, 'cm'))

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig3/fig3ab_legend_grp.pdf", 
    width = 12, height = 5, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(legend_grp))
dev.off() # Close the file

