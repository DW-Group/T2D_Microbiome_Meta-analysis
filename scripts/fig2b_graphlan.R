# prepare file for GraPhlAn
library(tidyverse)
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/analysis/abund0020/species_T2D_main.RData")
meta3extra <- meta3
meta_3extra <- meta_3
gdata::keep(meta3extra, meta_3extra, sure=T)

load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/analysis/abund0020/species_trend.RData")

taxa$feature <- taxa$X7
taxa$taxonomy <- gsub("k__", "", taxa$taxonomy)
taxa$taxonomy <- gsub("p__", "", taxa$taxonomy)
taxa$taxonomy <- gsub("c__", "", taxa$taxonomy)
taxa$taxonomy <- gsub("o__", "", taxa$taxonomy)
taxa$taxonomy <- gsub("f__", "", taxa$taxonomy)
taxa$taxonomy <- gsub("g__", "", taxa$taxonomy)
taxa$taxonomy <- gsub("s__", "", taxa$taxonomy)
taxa$taxonomy <- gsub("\\|", ".", taxa$taxonomy)
taxa <- taxa[(taxa$X7 %in% spe_final), ]

##################################################################
#         fixed-meta adjusting for age+sex+bmi+metf
##################################################################
extraspe <- setdiff(meta_3extra[meta_3extra$qval.fdr<0.1,]$feature,
                    meta_3[meta_3$qval.fdr<0.1,]$feature)

# replace those extra species with results from ConT2D model
pool0 <- rbind(meta_3[!(meta_3$feature %in% extraspe), ],
               meta_3extra[meta_3extra$feature %in% extraspe, ]) %>%
  mutate(lab=gsub("s__", "", feature)) %>%
  arrange(coef) %>%
  mutate(lab = factor(lab, levels = lab))
pool0 <- merge(taxa[, -8], pool0, by="feature") %>% filter(k>=4)

# tree file
write.table(pool0[,"taxonomy"],"/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/species_fig2b.txt",
            sep="", quote=F, row.names=FALSE, col.names=F)

# prepare annotation for rings
qin0 <- meta3$maaslin_fits$`Qin_2012 (CHN)` %>% mutate(study="Qin_2012 (CHN)", ring=8) %>% 
  filter(feature %in% pool0$feature)
karlsson0 <- meta3$maaslin_fits$`Karlsson_2013 (SWE)` %>% mutate(study="Karlsson_2013 (SWE)", ring=6) %>% 
  filter(feature %in% pool0$feature)
mbs0 <- meta3$maaslin_fits$`NHSII (USA)` %>% mutate(study="NHSII (USA)", ring=7) %>% 
  filter(feature %in% pool0$feature)
mlvs0 <- meta3$maaslin_fits$`HPFS (USA)` %>% mutate(study="HPFS (USA)", ring=5) %>% 
  filter(feature %in% pool0$feature)
zhong0 <- meta3$maaslin_fits$`Zhong_2019 (CHN)` %>% mutate(study="Zhong_2019 (CHN)", ring=12) %>% 
  filter(feature %in% pool0$feature)
direct0 <- meta3$maaslin_fits$`DIRECT-PLUS (ISR)` %>% mutate(study="DIRECT-PLUS (ISR)", ring=3) %>% 
  filter(feature %in% pool0$feature)
sol0 <- meta3$maaslin_fits$`SOL (USA)` %>% mutate(study="SOL (USA)", ring=9) %>% 
  filter(feature %in% pool0$feature)
wua0 <- meta3$maaslin_fits$`Wu_2020a (SWE)` %>% mutate(study="Wu_2020a (SWE)", ring=10) %>% 
  filter(feature %in% pool0$feature)
wub0 <- meta3$maaslin_fits$`Wu_2020b (SWE)` %>% mutate(study="Wu_2020b (SWE)", ring=11) %>% 
  filter(feature %in% pool0$feature)
fromentin0 <- meta3$maaslin_fits$`Fromentin_2022 (DEU_DNK_FRA)` %>% mutate(study="Fromentin_2022 (DEU/DNK/FRA)", ring=4) %>% 
  filter(feature %in% pool0$feature)

qin00 <- meta3extra$maaslin_fits$`Qin_2012 (CHN)` %>% mutate(study="Qin_2012 (CHN)", ring=8) %>% 
  filter(feature %in% pool0$feature)
karlsson00 <- meta3extra$maaslin_fits$`Karlsson_2013 (SWE)` %>% mutate(study="Karlsson_2013 (SWE)", ring=6) %>% 
  filter(feature %in% pool0$feature)
mbs00 <- meta3extra$maaslin_fits$`NHSII (USA)` %>% mutate(study="NHSII (USA)", ring=7) %>% 
  filter(feature %in% pool0$feature)
mlvs00 <- meta3extra$maaslin_fits$`HPFS (USA)` %>% mutate(study="HPFS (USA)", ring=5) %>% 
  filter(feature %in% pool0$feature)
zhong00 <- meta3extra$maaslin_fits$`Zhong_2019 (CHN)` %>% mutate(study="Zhong_2019 (CHN)", ring=12) %>% 
  filter(feature %in% pool0$feature)
direct00 <- meta3extra$maaslin_fits$`DIRECT-PLUS (ISR)` %>% mutate(study="DIRECT-PLUS (ISR)", ring=3) %>% 
  filter(feature %in% pool0$feature)
sol00 <- meta3extra$maaslin_fits$`SOL (USA)` %>% mutate(study="SOL (USA)", ring=9) %>% 
  filter(feature %in% pool0$feature)
wua00 <- meta3extra$maaslin_fits$`Wu_2020a (SWE)` %>% mutate(study="Wu_2020a (SWE)", ring=10) %>% 
  filter(feature %in% pool0$feature)
wub00 <- meta3extra$maaslin_fits$`Wu_2020b (SWE)` %>% mutate(study="Wu_2020b (SWE)", ring=11) %>% 
  filter(feature %in% pool0$feature)
fromentin00 <- meta3extra$maaslin_fits$`Fromentin_2022 (DEU_DNK_FRA)` %>% mutate(study="Fromentin_2022 (DEU/DNK/FRA)", ring=4) %>% 
  filter(feature %in% pool0$feature)

# replace those extra species with results from ConT2D model
qin0 <- rbind(qin0[!(qin0$feature %in% extraspe), ],
              qin00[qin00$feature %in% extraspe, ])
karlsson0 <- rbind(karlsson0[!(karlsson0$feature %in% extraspe), ],
                   karlsson00[karlsson00$feature %in% extraspe, ])
mbs0 <- rbind(mbs0[!(mbs0$feature %in% extraspe), ],
              mbs00[mbs00$feature %in% extraspe, ])
mlvs0 <- rbind(mlvs0[!(mlvs0$feature %in% extraspe), ],
               mlvs00[mlvs00$feature %in% extraspe, ])
zhong0 <- rbind(zhong0[!(zhong0$feature %in% extraspe), ],
                zhong00[zhong00$feature %in% extraspe, ])
direct0 <- rbind(direct0[!(direct0$feature %in% extraspe), ],
                 direct00[direct00$feature %in% extraspe, ])
sol0 <- rbind(sol0[!(sol0$feature %in% extraspe), ],
              sol00[sol00$feature %in% extraspe, ])
wua0 <- rbind(wua0[!(wua0$feature %in% extraspe), ],
              wua00[wua00$feature %in% extraspe, ])
wub0 <- rbind(wub0[!(wub0$feature %in% extraspe), ],
              wub00[wub00$feature %in% extraspe, ])
fromentin0 <- rbind(fromentin0[!(fromentin0$feature %in% extraspe), ],
                    fromentin00[fromentin00$feature %in% extraspe, ])

karlsson0 <- merge(taxa[, c("feature", "taxonomy")], karlsson0, by="feature")
qin0 <- merge(taxa[, c("feature", "taxonomy")], qin0, by="feature")
mbs0 <- merge(taxa[, c("feature", "taxonomy")], mbs0, by="feature")
mlvs0 <- merge(taxa[, c("feature", "taxonomy")], mlvs0, by="feature")
zhong0 <- merge(taxa[, c("feature", "taxonomy")], zhong0, by="feature")
direct0 <- merge(taxa[, c("feature", "taxonomy")], direct0, by="feature")
sol0 <- merge(taxa[, c("feature", "taxonomy")], sol0, by="feature")
wua0 <- merge(taxa[, c("feature", "taxonomy")], wua0, by="feature")
wub0 <- merge(taxa[, c("feature", "taxonomy")], wub0, by="feature")
fromentin0 <- merge(taxa[, c("feature", "taxonomy")], fromentin0, by="feature")
pool0$study="Pooled"
pool0$ring=2

rings0 <- rbind(pool0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")],
                qin0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")])
rings0 <- rbind(rings0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")],
                karlsson0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")])
rings0 <- rbind(rings0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")],
                mbs0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")])
rings0 <- rbind(rings0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")],
                mlvs0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")])
rings0 <- rbind(rings0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")],
                zhong0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")])
rings0 <- rbind(rings0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")],
                direct0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")])
rings0 <- rbind(rings0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")],
                sol0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")])
rings0 <- rbind(rings0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")],
                wua0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")])
rings0 <- rbind(rings0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")],
                wub0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")])
rings0 <- rbind(rings0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")],
                fromentin0[, c("taxonomy", "feature", "study", "coef", "stderr", "pval", "ring")])

rings0$feature <- gsub("s__", "", rings0$feature)
rings0$tvalue <- rings0$coef/rings0$stderr
summary(rings0$tvalue)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -4.11226 -0.81399 -0.01973 -0.00162  0.81968  4.91039     226 

rings0$color <- "ring_color"

library(dichromat)
colfunc_pos <- colorRampPalette(c("#bd0026", "white"))
colfunc_neg <- colorRampPalette(c("#08519c", "white"))
colfunc_pos(12)
# [1] "#BD0026" "#C31739" "#C92E4D" "#CF4561" "#D55C74" "#DB7388" "#E18B9C" "#E7A2B0" "#EDB9C3"
# [10] "#F3D0D7" "#F9E7EB" "#FFFFFF"
colfunc_neg(12)
# [1] "#08519C" "#1E60A5" "#3470AE" "#4B80B7" "#6190C0" "#78A0C9" "#8EAFD2" "#A5BFDB" "#BBCFE4"
# [10] "#D2DFED" "#E8EFF6" "#FFFFFF"

rings0 <- rings0 %>% mutate(group=case_when((tvalue<=(-4)) ~ "#3470AE",
                                            (tvalue>(-4) & tvalue<=(-3.5)) ~ "#4B80B7",
                                            (tvalue>(-3.5) & tvalue<=(-3)) ~ "#6190C0",
                                            (tvalue>(-3) & tvalue<=(-2.5)) ~ "#78A0C9",
                                            (tvalue>(-2.5) & tvalue<=(-2)) ~ "#8EAFD2",
                                            (tvalue>(-2) & tvalue<=(-1.5)) ~ "#A5BFDB",
                                            (tvalue>(-1.5) & tvalue<=(-1)) ~ "#BBCFE4",
                                            (tvalue>(-1) & tvalue<=(-0.5)) ~ "#D2DFED",
                                            (tvalue>(-0.5) & tvalue<=( 0)) ~ "#E8EFF6",
                                            (tvalue>( 0) & tvalue<=( 0.5)) ~ "#F9E7EB",
                                            (tvalue>( 0.5) & tvalue<=( 1)) ~ "#F3D0D7",
                                            (tvalue>( 1) & tvalue<=( 1.5)) ~ "#EDB9C3",
                                            (tvalue>( 1.5) & tvalue<=( 2)) ~ "#E7A2B0",
                                            (tvalue>( 2) & tvalue<=( 2.5)) ~ "#E18B9C",
                                            (tvalue>( 2.5) & tvalue<=( 3)) ~ "#DB7388",
                                            (tvalue>( 3) & tvalue<=( 3.5)) ~ "#D55C74",
                                            (tvalue>( 3.5) & tvalue<=( 4)) ~ "#CF4561",
                                            (tvalue>( 4) & tvalue<=( 4.5)) ~ "#C92E4D",
                                            (tvalue>( 4.5) & tvalue<=( 5)) ~ "#C31739"))

# significant species
sig_spe <- pool0 %>% filter(qval.fdr<0.1)
sig_spe$X1 <- gsub("k__", "", sig_spe$X1)
sig_spe$X2 <- gsub("p__", "", sig_spe$X2)
sig_spe$X3 <- gsub("c__", "", sig_spe$X3)
sig_spe$X4 <- gsub("o__", "", sig_spe$X4)
sig_spe$X5 <- gsub("f__", "", sig_spe$X5)
sig_spe$X6 <- gsub("g__", "", sig_spe$X6)
sig_spe$feature <- gsub("s__", "", sig_spe$feature)

# ring annotation file
rings0 <- rings0[(rings0$feature %in% sig_spe$feature), ]
write.table(rings0[!is.na(rings0$group),c("taxonomy", "color", "ring", "group")],
            "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/species_fig2b_ring3_12.txt",
            sep="\t", quote=F, row.names=FALSE, col.names=F)

# phylum color
taxa$phylum <- paste0(gsub("k__", "", taxa$X1), ".", gsub("p__", "", taxa$X2))
taxa <- taxa %>% filter(feature %in% pool0$feature) %>% 
  mutate(phylum_col=case_when(phylum=="Bacteria.Firmicutes" ~ "#282274",
                              phylum=="Bacteria.Actinobacteria" ~ "#1BAD4C",
                              phylum=="Bacteria.Bacteroidetes" ~ "#A7AC36",
                              phylum=="Bacteria.Proteobacteria" ~ "#1AAAAA",
                              phylum=="Bacteria.Verrucomicrobia" ~ "#DAB0D2",
                              phylum=="Archaea.Euryarchaeota" ~ "#4A6EB6",
                              phylum=="Bacteria.Fusobacteria" ~ "#F03B20",
                              phylum=="Bacteria.Lentisphaerae" ~ "#FDAE61"))
taxa$color <- "ring_color"
taxa$ring <- 1
write.table(taxa[,c("taxonomy", "color", "ring", "phylum_col")],
            "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/species_fig2b_ring1.txt",
            sep="\t", quote=F, row.names=FALSE, col.names=F)

taxa$clade_color <- "clade_marker_color"
write.table(taxa[!duplicated(taxa$X2),c("phylum", "clade_color", "phylum_col")],
            "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/species_fig2b_phy_color.txt",
            sep="\t", quote=F, row.names=FALSE, col.names=F)

# grey color to highlight pooled ring
pool0$color <- "ring_color"
pool0$grey <- "#c8c8c8"
pool0$grey2 <- "#888888"

write.table(pool0[pool0$qval.fdr>=0.1,c("taxonomy", "color", "ring", "grey2")],
            "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/species_fig2b_ring2.txt",
            sep="\t", quote=F, row.names=FALSE, col.names=F)

# species size & colour
sig_spe$clade_size <- "clade_marker_size"
sig_spe$size <- 75
write.table(sig_spe[,c("feature", "clade_size", "size")],
            "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/species_fig2b_spe_size.txt",
            sep="\t", quote=F, row.names=FALSE, col.names=F)
sig_spe$clade_color <- "clade_marker_color"
sig_spe <- sig_spe %>% mutate(color=case_when(X2=="Firmicutes" ~ "#282274",
                                              X2=="Actinobacteria" ~ "#1BAD4C",
                                              X2=="Bacteroidetes" ~ "#A7AC36",
                                              X2=="Proteobacteria" ~ "#1AAAAA",
                                              X2=="Verrucomicrobia" ~ "#DAB0D2",
                                              X2=="Euryarchaeota" ~ "#4A6EB6",
                                              X2=="Fusobacteria" ~ "#F03B20",
                                              X2=="Lentisphaerae" ~ "#FDAE61"))
write.table(sig_spe[,c("feature", "clade_color", "color")],
            "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/species_fig2b_spe_color.txt",
            sep="\t", quote=F, row.names=FALSE, col.names=F)

# genus annotation
sig_spe$annotation <- "annotation"
sig_spe$annotation_color <- "annotation_background_color"
sig_spe$annotation_font_size <- "annotation_font_size"
sig_spe$size <- 17

sig_spe$genus <- paste0(sig_spe$X1, ".",
                        sig_spe$X2, ".",
                        sig_spe$X3, ".",
                        sig_spe$X4, ".",
                        sig_spe$X5, ".",
                        sig_spe$X6)
write.table(sig_spe[!duplicated(sig_spe$X6),c("genus", "annotation", "X6")],
            "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/species_fig2b_genus_annot.txt",
            sep="\t", quote=F, row.names=FALSE, col.names=F)
write.table(sig_spe[!duplicated(sig_spe$X6),c("genus", "annotation_color", "color")],
            "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/species_fig2b_genus_color.txt",
            sep="\t", quote=F, row.names=FALSE, col.names=F)
write.table(sig_spe[!duplicated(sig_spe$X6),c("genus", "annotation_font_size", "size")],
            "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/species_fig2b_genus_size.txt",
            sep="\t", quote=F, row.names=FALSE, col.names=F)

# species abundance
sum(!is.na(species_T2D[1,])) # 4617
sum(species_T2D[1,], na.rm=T)
species_T2D2 <- species_T2D[rownames(species_T2D) %in% pool0$feature, ]
identical(pool0$feature, rownames(species_T2D2)) # TRUE
pool0$relab <- rowMeans(species_T2D2, na.rm=T)*100
summary(pool0$relab)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.003472 0.028925 0.101063 0.493818 0.403990 8.232267

pool0$relab_ring <- "ring_height"
pool0$ring2 <- "13"
write.table(pool0[,c("taxonomy", "relab_ring", "ring2", "relab")],
            "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/species_fig2b_relab.txt",
            sep="\t", quote=F, row.names=FALSE, col.names=F)
write.table(pool0[,c("taxonomy", "color", "ring2", "grey")],
            "/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/species_fig2b_relab_col.txt",
            sep="\t", quote=F, row.names=FALSE, col.names=F)

save.image("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/fig2b_graphlan.RData")

# spe_relab <- merge(spe_relab, spe_all[, c("species", "present")], by="species")
# spe_relab2 <- spe_relab %>% filter(species %in% meta_3$feature) %>% arrange(desc(relab))
# spe_relab2$rank <- 1:205
# spe_relab2$species <- gsub("s__", "", spe_relab2$species)
# spe_relab3 <- spe_relab2 %>% filter(species %in% gsub("s__", "", pool0$feature))
# View(spe_relab2[spe_relab2$species %in% sig_spe$feature,])
# length(intersect(spe_relab2[1:130,]$species, 
#                  sig_spe$feature))
# 
# length(intersect(spe_relab2[1:130,]$species,
#                  gsub("s__", "", pool0$feature)))
# 
# 
# sum(spe_relab2$present>=5) # 152
