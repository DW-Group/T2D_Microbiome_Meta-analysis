library(tidyverse)
library(MMUPHin)

setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/figS6")

load(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/analysis/abund0020/metf_species/metf_species.RData")

gdata::keep(meta_metf, meta_metf_results, spe_metf, sure=T)
load(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/analysis/abund0020/species_trend.RData")

spe_bmi <- meta_1[meta_1$qval.fdr<0.1, ]$feature # 46
spe_bmimetf <- meta_3[meta_3$qval.fdr<0.1, ]$feature # 14
spe_all <- union(spe_bmi, spe_bmimetf) # 46

meta_1_46 <- meta_1 %>% filter(feature %in% spe_all) %>% arrange(feature)
meta_3_46 <- meta_3 %>% filter(feature %in% spe_all) %>% arrange(feature)
identical(meta_1_46$feature, meta_3_46$feature) # TRUE

check_coef <- data.frame(species=spe_all,
                         coef_nometf=meta_1_46$coef,
                         p_nometf=meta_1_46$pval,
                         fdr_nometf=meta_1_46$qval.fdr,
                         coef_metf=meta_3_46$coef,
                         p_metf=meta_3_46$pval,
                         fdr_metf=meta_3_46$qval.fdr)
check_coef$per <- (check_coef$coef_metf-check_coef$coef_nometf)/check_coef$coef_nometf
check_coef$sign <- ifelse(check_coef$coef_metf * check_coef$coef_nometf>0, 1, 0)  

check_coef$metf10 <- ifelse(check_coef$species %in% 
                              (intersect(setdiff(spe_bmi, spe_bmimetf),
                                         meta_metf_results[meta_metf_results$fdr<0.1, ]$feature)), 
                            1, 0)
check_coef$metf4 <- ifelse(check_coef$species %in% 
                             (intersect(spe_bmimetf,
                                        meta_metf_results[meta_metf_results$fdr<0.1, ]$feature)), 
                           1, 0)

length(setdiff(spe_bmi, spe_bmimetf)) # 32
length(intersect(setdiff(spe_bmi, spe_bmimetf),
                 meta_metf_results[meta_metf_results$fdr<0.1, ]$feature)) # 10/32
# [1] "s__Acidaminococcus_intestini"       "s__Bilophila_wadsworthia"          
# [3] "s__Dorea_longicatena"               "s__Eubacterium_eligens"            
# [5] "s__Eubacterium_rectale"             "s__Faecalibacterium_prausnitzii"   
# [7] "s__Firmicutes_bacterium_CAG_95"     "s__Fusicatenibacter_saccharivorans"
# [9] "s__Intestinibacter_bartlettii"      "s__Roseburia_sp_CAG_309"

length(intersect(spe_bmimetf,
                 meta_metf_results[meta_metf_results$fdr<0.1, ]$feature)) # 2
# [1] "s__Clostridium_bolteae"    "s__Clostridium_sp_CAG_167" "s__Escherichia_coli"      
# [4] "s__Roseburia_sp_CAG_182" 

spe10 <- intersect(setdiff(spe_bmi, spe_bmimetf),
                   meta_metf_results[meta_metf_results$fdr<0.1, ]$feature)
spe4 <- intersect(spe_bmimetf,
                  meta_metf_results[meta_metf_results$fdr<0.1, ]$feature)

meta_1_10 <- meta_1[meta_1$feature %in% spe10, ]
meta_3_10 <- meta_3[meta_3$feature %in% spe10, ]
identical(meta_1_10$feature, meta_3_10$feature) # TRUE
sum(abs(meta_1_10$coef) > abs(meta_3_10$coef))
# 10 species, attenuated after adjusting for metformin

meta_1_4 <- meta_1[meta_1$feature %in% spe4, ]
meta_3_4 <- meta_3[meta_3$feature %in% spe4, ]
identical(meta_1_4$feature, meta_3_4$feature) # TRUE
sum(abs(meta_1_4$coef) > abs(meta_3_4$coef))
# 4 species, attenuated after adjusting for metformin

######################################################
#        10 species confounded by metf
######################################################
metf_plot10 <- meta_metf_results[meta_metf_results$feature %in% spe10,] %>% arrange(coef)

# data for p0_1
species0 <- metf_plot10 %>% mutate(lab = gsub("s__", "", feature))
species0$lab <- gsub("_", " ", species0$lab)
species0$lab <- factor(species0$lab, levels=species0$lab)
species0$lci <- species0$coef-species0$stderr
species0$uci <- species0$coef+species0$stderr

# data for p0_2
species_0 <- metf_plot10$feature
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
studies0 <- studies0 %>% mutate(star1=case_when(qval<0.05 ~ "#"),
                                star2=case_when(qval>=0.05 &  qval<0.1 ~ "*"))
studies0$feature <- factor(studies0$feature, levels=species0$feature)

# coef
summary(studies0$coef)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -5.2891 -0.7670 -0.3621 -0.2991  0.3853  2.2577 
sum(studies0$coef<(-3)) # 2
sum(studies0$coef>3) # 1
studies0$coef[studies0$coef>=3] <- 3
studies0$coef[studies0$coef<=(-3)] <- (-3)

summary(species0$lci)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -2.4189 -1.1317 -0.8372 -0.7623 -0.6354  1.1044
summary(species0$uci)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.2919 -0.6602 -0.4830 -0.2540 -0.3012  1.7961

species0_10 <- species0
studies0_10 <- studies0

library(ggplot2)
p10_1 <- ggplot(species0_10, aes(x=coef, y=lab, xmin=lci, xmax=uci)) +
  geom_errorbar(width=0.3, size=1) +
  geom_point(aes(color=fdr), size=3) +
  scale_color_gradient(low="#56B1F7", high="#132B43", 
                       limits=c(0,0.1), breaks=c(0, 0.1)) +
  geom_vline(xintercept=0, linetype="dashed", size=1.2) +
  scale_x_continuous(breaks=c(-2,-1,0,1,2), limits=c(-2.5, 2)) +
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
p10_2 <- ggplot(studies0_10, aes(x=study, y=feature)) +
  geom_point(aes(fill=coef), size=10, shape=21, stroke=0) +
  scale_fill_gradientn(
    colors=c("#08519c","white","#bd0026"),
    values=rescale(c(-3,0,3)),
    breaks=c(-3,0,3),
    limits=c(-3,3)) +
  geom_text(aes(label=star1), color="black", size=5, nudge_y = 0.0) +
  geom_text(aes(label=star2), color="black", size=7, nudge_y = -0.2) +
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
data1 <- meta_1[meta_1$feature %in% spe10, ]
data3 <- meta_3[meta_3$feature %in% spe10, ]
data1$model <- "Not adjusting for metf"
data3$model <- "Adjusting for metf"
heatmap <- rbind(data1[, c("feature", "model", "coef", "pval", "qval.fdr")],
                 data3[, c("feature", "model", "coef", "pval", "qval.fdr")])
heatmap <- heatmap %>% mutate(star1=case_when(qval.fdr<0.05 ~ "#"),
                              star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ ""))
heatmap$feature <- factor(heatmap$feature, levels=species0_10$feature)
summary(heatmap$coef)
#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# -0.565075 -0.422134 -0.020412  0.008045  0.406490  0.681908

sum(heatmap$coef>=0.5) # 1
sum(heatmap$coef<(-0.5)) # 1
heatmap$coef[heatmap$coef>=0.5] <- 0.5
heatmap$coef[heatmap$coef<(-0.5)] <- (-0.5)
heatmap$model <- factor(heatmap$model, levels=c("Not adjusting for metf",
                                                "Adjusting for metf"))

heatmap_10 <- heatmap
p10_3 <- ggplot(heatmap_10, aes(x=model, y=feature)) +
  geom_point(aes(color=coef), size=10, shape=15) +
  scale_color_gradient2(low="#08519c",
                        mid="white",
                        high="#bd0026",
                        midpoint = 0,
                        guide = "colourbar",
                        breaks=c(-0.5,0,0.5),
                        limits=c(-0.5,0.5)) +
  geom_text(aes(label=star1), color="black", size=5, nudge_y = 0.0) +
  geom_text(aes(label=star2), color="black", size=7, nudge_y = -0.2) +
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
#    4 species associated with both T2D and metf
######################################################
metf_plot4 <- meta_metf_results[meta_metf_results$feature %in% spe4,] %>% arrange(coef)

# data for p0_1
species0 <- metf_plot4 %>% mutate(lab = gsub("s__", "", feature))
species0$lab <- gsub("_", " ", species0$lab)
species0$lab <- factor(species0$lab, levels=species0$lab)
species0$lci <- species0$coef-species0$stderr
species0$uci <- species0$coef+species0$stderr

# data for p0_2
species_0 <- metf_plot4$feature
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
studies0 <- studies0 %>% mutate(star1=case_when(qval<0.05 ~ "#"),
                                star2=case_when(qval>=0.05 &  qval<0.1 ~ "*"))
studies0$feature <- factor(studies0$feature, levels=species0$feature)

# coef
summary(studies0$coef)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -2.50005 -0.80033 -0.06581  0.35740  1.35320  4.82567 
sum(studies0$coef<(-3)) # 0
sum(studies0$coef>3) # 2
studies0$coef[studies0$coef>=3] <- 3
studies0$coef[studies0$coef<=(-3)] <- (-3)

summary(species0$lci)

summary(species0$uci)

species0_4 <- species0
studies0_4 <- studies0
p4_1 <- ggplot(species0_4, aes(x=coef, y=lab, xmin=lci, xmax=uci)) +
  geom_errorbar(width=0.3, size=1) +
  geom_point(aes(color=fdr), size=3) +
  scale_color_gradient(low="#56B1F7", high="#132B43", 
                       limits=c(0,0.1), breaks=c(0, 0.1)) +
  geom_vline(xintercept=0, linetype="dashed", size=1.2) +
  scale_x_continuous(breaks=c(-1,0,1,2), limits=c(-1.2, 2.2)) +
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

p4_2 <- ggplot(studies0, aes(x=study, y=feature)) +
  geom_point(aes(fill=coef), size=10, shape=21, stroke=0) +
  scale_fill_gradientn(
    colors=c("#08519c","white","#bd0026"),
    values=rescale(c(-3,0,3)),
    breaks=c(-3,0,3),
    limits=c(-3,3)) +
  geom_text(aes(label=star1), color="black", size=5, nudge_y = 0.0) +
  geom_text(aes(label=star2), color="black", size=9, nudge_y = -0.2) +
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
data1 <- meta_1[meta_1$feature %in% spe4, ]
data3 <- meta_3[meta_3$feature %in% spe4, ]
data1$model <- "Not adjusting for metf"
data3$model <- "Adjusting for metf"
heatmap <- rbind(data1[, c("feature", "model", "coef", "pval", "qval.fdr")],
                 data3[, c("feature", "model", "coef", "pval", "qval.fdr")])
heatmap <- heatmap %>% mutate(star1=case_when(qval.fdr<0.05 ~ "#"),
                              star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "*"))
heatmap$feature <- factor(heatmap$feature, levels=species0_4$feature)
summary(heatmap$coef)

sum(heatmap$coef>=0.5) # 1
sum(heatmap$coef<(-0.5)) # 1
heatmap$coef[heatmap$coef>=0.5] <- 0.5
heatmap$coef[heatmap$coef<(-0.5)] <- (-0.5)
heatmap$model <- factor(heatmap$model, levels=c("Not adjusting for metf",
                                                "Adjusting for metf"))

heatmap_4 <- heatmap
p4_3 <- ggplot(heatmap_4, aes(x=model, y=feature)) +
  geom_point(aes(color=coef), size=10, shape=15) +
  scale_color_gradient2(low="#08519c",
                        mid="white",
                        high="#bd0026",
                        midpoint = 0,
                        guide = "colourbar",
                        breaks=c(-0.5,0,0.5),
                        limits=c(-0.5,0.5)) +
  geom_text(aes(label=star1), color="black", size=5, nudge_y = 0.0) +
  geom_text(aes(label=star2), color="black", size=9, nudge_y = -0.2) +
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

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/figS6/fig_s6b.pdf", 
    width = 11.9, height = 9.55, onefile = F) # Open a new pdf file
egg::ggarrange(p10_1, p10_2, p10_3, 
               p4_1, p4_2, p4_3, 
               nrow=2, ncol=3, 
               widths = c(2,1.6,0.48),
               heights = c(11,4.48))
dev.off() # Close the file
