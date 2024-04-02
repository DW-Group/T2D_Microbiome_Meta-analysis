setwd("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2")

library(tidyverse)
load("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/analysis/abund0020/species_trend.RData")
meta_1trend <- meta_1
meta_2trend <- meta_2
meta_3trend <- meta_3

gdata::keep(meta_1trend, meta_2trend, meta_3trend, taxa, sure=T)

spe_all <- Reduce(union, list(
  meta_1trend[meta_1trend$qval.fdr<0.1,]$feature,
  meta_2trend[meta_2trend$qval.fdr<0.1,]$feature,
  meta_3trend[meta_3trend$qval.fdr<0.1,]$feature
))

data1trend <- meta_1trend[meta_1trend$feature %in% spe_all, ]
data2trend <- meta_2trend[meta_2trend$feature %in% spe_all, ]
data3trend <- meta_3trend[meta_3trend$feature %in% spe_all, ]
data1trend$model <- "age + sex + BMI"
data2trend$model <- "age + sex + metf"
data3trend$model <- "age + sex + BMI + metf"

heatmap_trend <- rbind(
  data1trend[, c("feature", "model", "coef", "pval", "qval.fdr")],
  data2trend[, c("feature", "model", "coef", "pval", "qval.fdr")],
  data3trend[, c("feature", "model", "coef", "pval", "qval.fdr")]
)

heatmap_trend <- heatmap_trend %>% mutate(star1=case_when(qval.fdr<0.05 ~ "#"),
                                          star2=case_when(qval.fdr>=0.05 &  qval.fdr<0.1 ~ "*"))
heatmap_trend$species <- gsub("s__", "", heatmap_trend$feature)
heatmap_trend$species <- gsub("_", " ", heatmap_trend$species)

colnames(taxa)[8] <- "feature"
heatmap_trend <- merge(taxa[, c("feature", "X2")], heatmap_trend, by="feature", all.y=T)
heatmap_trend$X2 <- gsub("p__", "", heatmap_trend$X2)

summary(heatmap_trend$coef)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.56507 -0.22588 -0.12095 -0.01969  0.20037  0.68191 

sum(heatmap_trend$coef>0.7) # 0
sum(heatmap_trend$coef<(-0.5)) # 3
heatmap_trend$coef[heatmap_trend$coef>=0.7] <- 0.7
heatmap_trend$coef[heatmap_trend$coef<(-0.5)] <- (-0.5)

sum(data1trend$coef>0) # 24
data1trend <- merge(taxa[, c("feature", "X2")], data1trend, by="feature", all.y=T)
data1trend$X2 <- gsub("p__", "", data1trend$X2)
data1trend$species <- gsub("s__","",data1trend$feature)
data1trend$species <- gsub("_"," ",data1trend$species)

data1trend$sign <- ifelse(data1trend$coef>0,1,0)
data1trend <- data1trend %>% arrange(sign, desc(X2), desc(species))

heatmap_trend$sign <- ifelse((heatmap_trend$feature %in% data1trend[data1trend$coef>0,]$feature),1,0)
heatmap_trend$model <- factor(heatmap_trend$model, levels=c("age + sex + BMI", 
                                                            "age + sex + metf", 
                                                            "age + sex + BMI + metf"))

heatmap_trend1 <- heatmap_trend %>% filter(sign==1)
heatmap_trend1$species <- factor(heatmap_trend1$species, levels=data1trend[data1trend$coef>0,]$species)
heatmap_trend2 <- heatmap_trend %>% filter(sign==0)
heatmap_trend2$species <- factor(heatmap_trend2$species, levels=data1trend[data1trend$coef<0,]$species)

table(heatmap_trend[heatmap_trend$qval.fdr<0.1,]$model)
# age + sex + BMI       age + sex + metf age + sex + BMI + metf 
# 46                                  40                     14 

beta <- data.frame(mod1trend=meta_1trend$coef,
                   mod2trend=meta_2trend$coef,
                   mod3trend=meta_3trend$coef)
cor(beta)
#           mod1trend mod2trend mod3trend
# mod1trend 1.0000000 0.8981223 0.9071312
# mod2trend 0.8981223 1.0000000 0.9690377
# mod3trend 0.9071312 0.9690377 1.0000000

phycol <- c("#1BAD4C", "#A7AC36", "#4A6EB6", "#282274", "#F1B16E", "#1AAAAA")
names(phycol) <- c("Actinobacteria", "Bacteroidetes", "Euryarchaeota", "Firmicutes",
                   "Lentisphaerae", "Proteobacteria")

library(ggplot2)
heat1_0 <- ggplot(heatmap_trend1, aes(x=1, y=species)) +
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

heat1_trend <- ggplot(heatmap_trend1, aes(x=model, y=species)) +
  geom_tile(aes(fill=coef)) +
  scale_fill_gradient2(low="#08519c",
                       mid="white",
                       high="#bd0026",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.5,0,0.7),
                       limits=c(-0.5,0.7)) +
  scale_x_discrete(position = "top") +
  geom_text(aes(label=star1), color="black", size=5, nudge_y = 0.05) +
  geom_text(aes(label=star2), color="black", size=9, nudge_y = -0.3) +
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

heat2_0 <- ggplot(heatmap_trend2, aes(x=1, y=species)) +
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

heat2_trend <- ggplot(heatmap_trend2, aes(x=model, y=species)) +
  geom_tile(aes(fill=coef)) +
  scale_fill_gradient2(low="#08519c",
                       mid="white",
                       high="#bd0026",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.5,0,0.7),
                       limits=c(-0.5,0.7)) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  geom_text(aes(label=star1), color="black", size=5, nudge_y = 0.05) +
  geom_text(aes(label=star2), color="black", size=9, nudge_y = -0.3) +
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
        legend.position = "none",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14),
        # legend.justification = c(2, 0), 
        legend.box.margin=margin(c(0,-100,0,0))) +
  guides(fill = guide_colourbar(title.position="top",
                                title.hjust = 0.5,
                                barheight = 1.5,
                                barwidth = 10))

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/fig2a_v2.pdf",
    width = 7.45, height = 24, onefile = F) # Open a new pdf file
egg::ggarrange(heat1_trend, heat1_0, 
               heat2_trend, heat2_0, nrow=2, ncol=2,
               widths = c(1,0.1), heights=c(24,34))
dev.off() # Close the file

# legend
trendlegend <- ggplot(heatmap_trend2, aes(x=model, y=species)) +
  geom_tile(aes(fill=coef)) +
  scale_fill_gradient2(low="#08519c",
                       mid="white",
                       high="#bd0026",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.5,0,0.7),
                       limits=c(-0.5,0.7)) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "right") +
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
        legend.position = "top",
        legend.text = element_text(colour = "black", size = 14),
        legend.title = element_text(colour = "black", size = 14)) +
  guides(fill = guide_colourbar(title.position="top",
                                title.hjust = 0.5,
                                barheight = 1.5,
                                barwidth = 10))

pdf("/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/fig2a_legend.pdf",
    width = 3, height = 2, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(trendlegend))
dev.off() # Close the file

save.image(file="/n/holylabs/LABS/dwang_lab/Lab/T2D_meta/3.analysis/R2/fig2/fig2a.RData")

# mean/se for trend
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -3.4213 -2.4149  2.3084  0.6684  3.0459  5.1359 

# mean/se for T2DCon
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -3.6626 -2.4319  2.3659  0.4945  2.9224  4.3925
