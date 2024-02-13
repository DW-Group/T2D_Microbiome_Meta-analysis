setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R1/figS5")

library(tidyverse)
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/0323/species_trend.RData")
gdata::keep(meta_0, meta_1, meta_2, meta_3, meta_4, taxa, sure=T)

spe_all <- Reduce(union, list(
  meta_0[meta_0$qval.fdr<0.25,]$feature,
  meta_1[meta_1$qval.fdr<0.25,]$feature,
  meta_2[meta_2$qval.fdr<0.25,]$feature,
  meta_3[meta_3$qval.fdr<0.25,]$feature,
  meta_4[meta_4$qval.fdr<0.25,]$feature
))

data0 <- meta_0[meta_0$feature %in% spe_all, ]
data1 <- meta_1[meta_1$feature %in% spe_all, ]
data2 <- meta_2[meta_2$feature %in% spe_all, ]
data3 <- meta_3[meta_3$feature %in% spe_all, ]
data4 <- meta_4[meta_4$feature %in% spe_all, ]
data0$model <- "age + sex"
data1$model <- "age + sex + BMI"
data2$model <- "age + sex + metf"
data3$model <- "age + sex + BMI + metf"
data4$model <- "age + sex + BMI + metf + insul"

heatmap <- rbind(
  data0[, c("feature", "model", "coef", "pval", "qval.fdr")],
  data1[, c("feature", "model", "coef", "pval", "qval.fdr")],
  data2[, c("feature", "model", "coef", "pval", "qval.fdr")],
  data3[, c("feature", "model", "coef", "pval", "qval.fdr")],
  data4[, c("feature", "model", "coef", "pval", "qval.fdr")]
)
heatmap <- heatmap %>% mutate(star1=case_when(qval.fdr<0.05 ~ "+"),
                              star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "#"),
                              star3=case_when(qval.fdr>=0.1 &  qval.fdr<0.25 ~ "*"))
heatmap$species <- gsub("s__", "", heatmap$feature)
heatmap$species <- gsub("_", " ", heatmap$species)

colnames(taxa)[8] <- "feature"
heatmap <- merge(taxa[, c("feature", "X2")], heatmap, by="feature", all.y=T)
heatmap$X2 <- gsub("p__", "", heatmap$X2)

summary(heatmap$coef)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.660285 -0.164046  0.003085  0.026354  0.201126  0.717679

sum(heatmap$coef>0.7) # 2
heatmap$coef[heatmap$coef>0.7] <- 0.7

sum(data1$coef>0) # 52
data1 <- merge(taxa[, c("feature", "X2")], data1, by="feature", all.y=T)
data1$X2 <- gsub("p__", "", data1$X2)
data1$species <- gsub("s__","",data1$feature)
data1$species <- gsub("_"," ",data1$species)

data1 <- data1 %>% arrange(desc(X2), desc(species))

heatmap$sign <- ifelse((heatmap$feature %in% data1[data1$coef>0,]$feature),1,0)
heatmap$model <- factor(heatmap$model, levels=c(
  "age + sex",
  "age + sex + BMI",
  "age + sex + metf",
  "age + sex + BMI + metf",
  "age + sex + BMI + metf + insul"
))

heatmap1 <- heatmap %>% filter(sign==1)
heatmap1$species <- factor(heatmap1$species, levels=data1[data1$coef>0,]$species)
heatmap2 <- heatmap %>% filter(sign==0)
heatmap2$species <- factor(heatmap2$species, levels=data1[data1$coef<0,]$species)

table(heatmap[heatmap$qval.fdr<0.25,]$model)
# age + sex                age + sex + BMI               age + sex + metf 
# 101                             90                             63 
# age + sex + BMI + metf age + sex + BMI + metf + insul 
# 37                             30 

beta <- data.frame(mod0=meta_0$coef,
                   mod1=meta_1$coef,
                   mod2=meta_2$coef,
                   mod3=meta_3$coef,
                   mod4=meta_4$coef)
cor(beta)

phycol <- c("#1BAD4C", "#A7AC36", "#4A6EB6", "#282274", "#FDAE61", "#1AAAAA")
names(phycol) <- c("Actinobacteria", "Bacteroidetes", "Euryarchaeota", "Firmicutes",
                   "Lentisphaerae", "Proteobacteria")

library(ggplot2)
heat1_0 <- ggplot(heatmap1, aes(x=1, y=species)) +
  geom_tile(aes(fill=X2)) +
  labs(fill="Phylum") +
  scale_fill_manual(values=phycol) +
  xlab("") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "italic", color = "black", size = 14),
        axis.text.x = element_blank(),
        axis.title =element_blank(),
        legend.position = "none",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14)
  )

heat1 <- ggplot(heatmap1, aes(x=model, y=species)) +
  geom_tile(aes(fill=coef)) +
  scale_fill_gradient2(low="#08519c",
                       mid="white",
                       high="#bd0026",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.7,0,0.7),
                       limits=c(-0.7,0.7)) +
  scale_x_discrete(position = "top") +
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
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 14, colour = "black"),
        # axis.ticks.length=unit(0.07,"inch"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5)) +
  theme(legend.position = "none")

heat2_0 <- ggplot(heatmap2, aes(x=1, y=species)) +
  geom_tile(aes(fill=X2)) +
  scale_fill_manual(values=phycol) +
  scale_y_discrete(position = "right") +
  xlab("") +
  theme_bw() +
  theme(panel.border = element_rect(size=1.2), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "italic", color = "black", size = 14),
        axis.text.x = element_blank(),
        axis.title =element_blank(),
        legend.position = "none")

heat2 <- ggplot(heatmap2, aes(x=model, y=species)) +
  geom_tile(aes(fill=coef)) +
  scale_fill_gradient2(low="#08519c",
                       mid="white",
                       high="#bd0026",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.7,0,0.7),
                       limits=c(-0.7,0.7)) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  geom_text(aes(label=star1), color="black", size=5, nudge_y = 0.05) +
  geom_text(aes(label=star2), color="black", size=5, nudge_y = 0.0) +
  geom_text(aes(label=star3), color="black", size=7, nudge_y = -0.2) +
  labs(fill=expression(beta * " coefficient")) +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_line(colour = "black", size = 0.8),
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(angle = 45, vjust = 0, hjust = 0, size = 14, colour = "black"),
        axis.ticks.length=unit(0.04,"inch"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),
        legend.position = "bottom",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14),
        # legend.justification = c(2, 0), 
        legend.box.margin=margin(c(0,-100,0,0))) +
  guides(fill = guide_colourbar(title.position="top",
                                title.hjust = 0.5,
                                barheight = 1.5,
                                barwidth = 10))

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R1/figS5/figure_s5a.pdf",
    width = 7.7, height = 36.5, onefile = F) # Open a new pdf file
egg::ggarrange(heat1, heat1_0, heat2, heat2_0, nrow=2, ncol=2,
               widths = c(1, 0.1), heights=c(52,56))
dev.off() # Close the file

 
# # legend
# phylum <- ggplot(heatmap1, aes(x=1, y=species)) +
#   geom_tile(aes(fill=X2)) +
#   labs(fill="Phylum") +
#   scale_fill_manual(values=phycol) +
#   xlab("") +
#   theme_bw() +
#   theme(panel.border = element_rect(size=1.2), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), 
#         axis.ticks.y = element_line(colour = "black", size = 0.8),
#         axis.ticks.x = element_blank(),
#         axis.text.y = element_text(face = "italic", color = "black", size = 14),
#         axis.text.x = element_blank(),
#         axis.title =element_blank(),
#         legend.position = "bottom",
#         legend.text = element_text(colour = "black", size = 14),
#         legend.title = element_text(colour = "black", size = 14),
#         legend.spacing.y = unit(0.1, 'cm'),
#         legend.spacing.x = unit(0.15, 'cm')) +
#   guides(fill = guide_legend(title.position="top", nrow=3, byrow=TRUE))
# 
# pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s4/figure_s4a_legend.pdf",
#     width = 6, height = 3, onefile = F) # Open a new pdf file
# grid::grid.draw(ggpubr::get_legend(phylum))
# dev.off() # Close the file
# 
# # scatter plot
# p1 <- ggplot(aes(x=mod0, y=mod1), data=beta) +
#   geom_point() +
#   geom_hline(yintercept = 0, col="red") +
#   geom_vline(xintercept = 0, col="red") +
#   scale_x_continuous(breaks=c(-0.6,-0.3,0,0.3,0.6), limits=c(-0.68,0.72)) +
#   scale_y_continuous(breaks=c(-0.3,0,0.3,0.6), limits=c(-0.5, 0.72)) +
#   labs(x="Adjustment for\nage + sex",
#        y="Adjustment for\nage + sex + BMI",
#        title="") +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid = element_blank(),
#         axis.text = element_text(colour = "black", size=12),
#         axis.title = element_text(colour = "black", size=12),
#         panel.border = element_rect(colour = "black", size=1.2)) +
#   annotate(geom="text", x=-0.45, y=0.63, col="black", label="r=0.98", size=4.2)
# 
# p3 <- ggplot(aes(x=mod1, y=mod3), data=beta) +
#   geom_point() +
#   geom_hline(yintercept = 0, col="red") +
#   geom_vline(xintercept = 0, col="red") +
#   scale_x_continuous(breaks=c(-0.3,0,0.3,0.6), limits=c(-0.5, 0.72)) +
#   scale_y_continuous(breaks=c(-0.4,-0.2,0,0.2,0.4), limits=c(-0.4,0.54)) +
#   labs(x="Adjustment for\nage + sex + BMI",
#        y="Adjustment for\nage + sex + BMI + metf",
#        title="") +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid = element_blank(),
#         axis.text = element_text(colour = "black", size=12),
#         axis.title = element_text(colour = "black", size=12),
#         panel.border = element_rect(colour = "black", size=1.2)) +
#   annotate(geom="text", x=-0.3, y=0.45, col="black", label="r=0.91", size=4.2)
# 
# pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s4/figure_s4bc.pdf",
#     width = 6, height = 3, onefile = F) # Open a new pdf file
# egg::ggarrange(p1, p3, nrow=1, ncol=2)
# dev.off() # Close the file
# 
# save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/dropbox/fig_s4/figure_s4.RData")
# 
